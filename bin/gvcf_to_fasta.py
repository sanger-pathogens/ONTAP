#!/usr/bin/env python
from pysam import VariantFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import argparse, sys, os
import pandas as pd

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print(f"{message}\n\n")
        self.print_help()
        sys.exit(2)

    def is_valid_file(self, parser, arg):
        if not os.path.isfile(arg):
            parser.error(f"The file {arg} does not exist!")
        else:
            return arg


def argparser():
    description = """
    Tool to plot variants onto a reference region as described by a bedfile
    """
    parser = ParserWithErrors(description = description)
    # input files
    parser.add_argument("-r", "--reference_file", required=True,
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-v", "--gvcf_file", required=True,
                        help="input gVCF file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-b", "--bed_file", type=lambda x: parser.is_valid_file(parser, x), required=True,
                        help="BED file (TSV) defining regions (<name>\t<start>\t<end>)" )
    # general output
    parser.add_argument("-o", "--output_fasta_file_prefix", required=True,
                    help="prefix for output fasta file")
    parser.add_argument("-i", "--fasta_id",
                    default="auto", help="fasta header ID")
    # output file type
    parser.add_argument("-m", "--multifasta",
                    action="store_true", help="output one multifasta file per sample containing a sequence record per locus")
    parser.add_argument("-s", "--singlefasta",
                    action="store_true", help="output a fasta file per locus and per sample containing a single sequence record each")
    parser.add_argument("-ml", "--multi_locus",
                    action="store_true", help="output a multi-locus concatenated fasta containing a sequence record per sample")
    parser.add_argument("-w", "--whole_genome_fasta",
                    action="store_true", help="output a (multi)fasta file with whole-geneome sequence representation, one file per sample containing a sequence record per reference chromosome")
    # processing options
    parser.add_argument("-n", "--unknown_as_n",
                    action="store_true", help="represent genoytpe calls that are not supported with Ns; default is to use the reference sequence allele in that position")
    parser.add_argument("-g", "--gap_character", default= 'N', 
                        help="character to use for between amplicon gaps default - leave blank to only include variants (bit useless)")
    parser.add_argument("-q", "--min_alt_gt_qual", default= '1', type=int,
                        help="minimum genotype quality (GQ) for reporting an ALT variant in consensus sequence; recommended value: 5")
    parser.add_argument("-Q", "--min_ref_gt_qual", default= '1', type=int,
                        help="minimum genotype quality (GQ) for reporting a REF variant in consensus sequence; recommended value: 4")
    return parser

def get_variant_info(gvcf_file, min_ref_gt_qual, min_alt_gt_qual):
    """
    extracts variants from gVCF file

    gVCF file (file): gVCF file to extract variants from

    returns
    returns a ordereddict of variants from the input gVCF
    """
    variant_info = OrderedDict()
    with VariantFile(gvcf_file) as gvcf:
        for record in gvcf:
            position = record.pos - 1 # make it 0-based
            ref_allele = record.ref
            gq = record.samples.items()[0][1]['GQ']
            end_block = None
            chrom = record.chrom

            if record.alts is None:
                called_allele = ref_allele
            elif len(record.alts) == 1:
                if record.alts[0] == '<NON_REF>': # non intuitive from the allele being called NON_REF, but all genotype calls have this
                    if gq >= min_ref_gt_qual:
                        called_allele = ref_allele
                        end_block = record.stop
                    else:
                        continue
                else:
                    raise ValueError(f"unexpected single allele value for ALT field: {record.alts}")
            elif len(record.alts) == 2:
                if record.alts[1] != '<NON_REF>':
                    raise ValueError(f"unexpected multiple allele value for ALT field: {record.alts}")
                else:
                    if gq >= min_alt_gt_qual:
                        called_allele = record.alts[0] # actual alt non-ref variant
                    else:
                        continue
            else:
                raise ValueError(f"unexpected multiple allele value for ALT field: {record.alts}")
            variant_info[position] = (ref_allele, called_allele, end_block, chrom)
    return variant_info

def expand_variant_blocks(variants, ref_dict):
    expanded_variants = {}
    for key, value in variants.items():
        pos = key
        ref_allele, called_allele, end_block, chrom = value
        expanded_variants[pos] = value[:2]
        if not end_block is None:
            # detects definition of a genotype block
            contig_len = len(ref_dict[chrom])
            while pos < end_block:
                pos += 1
                if pos >= contig_len:   # guard: VCF end_block is 1-based, ref_dict is 0-based
                    break
                ref_allele = ref_dict[chrom][pos]
                expanded_variants[pos] = (ref_allele, ref_allele)
    
    return expanded_variants


def variants_in_range(bed_range, variants):
    """
    script to check if any variants fall in the bed region

    bed range (range): range of loci from the bed file
    variants (dict): base change and location of variant

    returns
    variants list (list): list of variants in the bed region
    """
    variants_list = []
    for key, value in variants.items():
        if key in bed_range:
            variants_list.append((key, value))
    return variants_list

def change_base_with_checks(sequence, position, base_change):
    """
    Change the base at a specified position in the sequence to a new base.

    Parameters:
    sequence (str): The biological sequence (e.g., DNA or RNA).
    position (int): The position (0-based index) where the base should be changed.
    base_change (tuple): Tuple containing old base and new one to replace it with.

    Returns:
    str: The modified sequence with the base changed.
    """
    old_base, new_base = base_change
    if position < 0 or position >= len(sequence):
        raise ValueError(f"Position {position} is out of range: {len(sequence)}")
    if old_base != sequence[position].upper() and sequence[position].upper() != "N":
        raise ValueError(f"base at index {position} was expected to be {old_base} in {sequence}")
    if new_base not in 'ACGTUNacgtun':
        print(position, base_change)
        raise ValueError(f"Variant {new_base} not an accepted base.")
    
    
    # Change the base at the specified position
    new_sequence = sequence[:position] + new_base + sequence[position+1:]
    return new_sequence

def calculate_gaps_to_add(gap_start_position, gap_end_position, gapcharacter):
    """
    takes a start and end position and creates a gap that length with a character of your choosing

    returns
    a list of gap characters n long depending on input
    """
    return [gapcharacter] * (gap_end_position - gap_start_position)

def extract_sequences_from_bed_and_include_called_genotypes(ref_dict, bed_df, variant_info, unknown_as_n, gap_character):
    """
    Generate background reference (or gaps) for loci in a bed file replacing reference with variants where they are found

    Parameters:
    ref_dict (dict of SeqRecords): The Reference biological sequence (e.g., DNA or RNA).
    bed_df (data.frame): The positions of the loci's for the above reference
    variant_info (file): A VCF file containing variant calls for your chosen sample against the reference
    unknown_as_n (str): A binary yes no if the reference should be replaced with unknown base characters
    gap_character (str): Where reads cannot support calling a genotype, a character to replace the reference base with (usually N)

    Returns:
    list: A list of Seq objects that are reflective of the loci from the reference with variant bases overwritten with their variants as called in the gVCF file
    """
    # Extract sequences
    extracted_sequences = []
    for index, row in bed_df.iterrows():
        chromosome = row['chromosome'] #not needed for this but good to ensure we are looking in the correct place? Multibed/multifasta ref later?
        start = row['start']
        end = row['end']

        if unknown_as_n:
            gaps = calculate_gaps_to_add(start, end, gap_character)
            sequence = Seq("".join(gaps))
        else:
            sequence = ref_dict[chromosome].seq[start:end].lower()
        
        variants = variants_in_range(range(start, end), variant_info)
        
        for variant in variants:
            adjusted_start = variant[0] - start
            sequence = change_base_with_checks(sequence, adjusted_start, variant[1])

        sequence_record = SeqRecord(sequence, id=f"{chromosome}_{start}_{end}", description=f"{start}_{end}")
        extracted_sequences.append(sequence_record)
    
    return extracted_sequences

def write_sequence(fasta_prefix, fasta_id, sequence_list, multi_locus=True, multifasta=False, singlefasta=False):
    """
    Write sequences to file

    fasta_prefix (path): prefix of fasta to write to
    multi_locus (bool): write multilocus concatenate output or not
    multifasta (bool): write multifasta output or not
    fasta id (str): ID for the fasta header
    sequence list (list): A list of SEQ records to be written to file either as a single fasta or multifasta

    TO DO : add whole-genome consensus sequence output

    returns
    writes files returns nothing
    """
    if multifasta:
        with open(f"{fasta_prefix}.fasta", 'w') as output:
            for i, sequence in enumerate(sequence_list, start=1):
                record = SeqRecord(Seq(sequence.seq), id = f"{fasta_id}_{i}", description = '')
                SeqIO.write(record, output, "fasta")
    if multi_locus:
        sequences = [str(sequence.seq) for sequence in sequence_list]
        master_record = "".join(sequences)
        with open(f"{fasta_prefix}_multi_locus.fasta", 'w') as output:
            record = SeqRecord(Seq(master_record), id = f"{fasta_id}_multi_locus", description = 'joined ref bed file regions with variants') 
            SeqIO.write(record, output, "fasta")        
    
    if singlefasta:
        for sequence in sequence_list:
            start, end = sequence.description.split("_")
            final_name = f"{fasta_prefix}_{start}_{end}"
            with open(f"{final_name}.fa", 'w') as output:
                    record = SeqRecord(Seq(sequence.seq), id = f"{fasta_id}_{sequence.id}", description = '')
                    SeqIO.write(record, output, "fasta")    

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

    #Reference genome as seqrecord
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.reference_file, "fasta"))

    # Read the bed file in
    locus_bed_df = pd.read_csv(args.bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end'], usecols=[0, 1, 2])

    called_genotypes = get_variant_info(args.gvcf_file, args.min_ref_gt_qual, args.min_alt_gt_qual)

    all_called_genotypes = expand_variant_blocks(called_genotypes, ref_dict)

    extracted_consensus_sequences = extract_sequences_from_bed_and_include_called_genotypes(ref_dict, locus_bed_df, all_called_genotypes, args.unknown_as_n, args.gap_character)

    if args.fasta_id == "auto":
        fasta_id = os.path.basename(args.gvcf_file).split('.')[0]
    else:
        fasta_id = args.fasta_id

    write_sequence(fasta_prefix=args.output_fasta_file_prefix, 
                   multi_locus=args.multi_locus, 
                   multifasta=args.multifasta, 
                   singlefasta=args.singlefasta, 
                   fasta_id=fasta_id, 
                   sequence_list=extracted_consensus_sequences)
    
    if args.whole_genome_fasta:
        wgs_bed_df = pd.DataFrame([(chrom, 0, len(seqrec)) for chrom, seqrec in ref_dict.items()], columns=locus_bed_df.columns)
        wg_consensus_sequences = extract_sequences_from_bed_and_include_called_genotypes(ref_dict, wgs_bed_df, all_called_genotypes, args.unknown_as_n, args.gap_character)
        write_sequence(fasta_prefix=f"{args.output_fasta_file_prefix}_wg",
                   multi_locus=False,
                   multifasta=args.multifasta, 
                   singlefasta=False, 
                   fasta_id=fasta_id, 
                   sequence_list=wg_consensus_sequences)