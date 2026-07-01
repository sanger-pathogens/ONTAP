#!/usr/bin/env nextflow

include { INDEX_REF } from './modules/samtools.nf'
include { MULTIQC } from './modules/multiqc.nf'

include { BASECALLING } from './subworkflows/basecalling.nf'
include {
    PRE_MAP_QC as PRE_MAP_QC_PRE_TRIM;
    PRE_MAP_QC as PRE_MAP_QC_POST_TRIM;
} from './subworkflows/pre_map_qc.nf'
include { PROCESS_FILTER_READS } from './subworkflows/process_filter_reads.nf'
include { MAPPING } from './subworkflows/mapping.nf'
include { FILTER_BAM } from './subworkflows/post_map_filtering.nf'
include {
    POST_MAP_QC;
    POST_FILTER_QC;
} from './subworkflows/post_map_qc.nf'
include { CALL_VARIANTS } from './subworkflows/variant_calling.nf'

process CAT_FASTQ {
    tag "${meta.ID}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.ID}_merged.fastq.gz")

    script:
    """
    cat ${reads.sort().join(' ')} > ${meta.ID}_merged.fastq.gz
    """
}

workflow {
    // Safe defaults (also useful if config missing values)
    params.outdir            = params.outdir ?: "results"
    params.save_fastqs       = params.save_fastqs ?: true
    params.keep_sorted_bam   = params.keep_sorted_bam ?: true
    params.keep_bam_files    = params.keep_bam_files ?: false
    params.monochrome_logs   = params.monochrome_logs ?: false
    params.dorado_local_path = params.dorado_local_path ?: ""
    def basecall_raw = params.basecall
    def do_basecall = (basecall_raw != null) && basecall_raw.toString().trim().equalsIgnoreCase('true')

    def barcode_kit = params.barcode_kit ?: (
        (params.barcode_kit_name instanceof List && params.barcode_kit_name.size() > 0)
            ? params.barcode_kit_name[0]
            : params.barcode_kit_name
    )

    if (do_basecall) {
        if (!params.raw_read_dir) exit 1, "basecall=true but --raw_read_dir is not set"
        if (!params.basecall_model) exit 1, "basecall=true but --basecall_model is not set"
    } else {
        if (!params.fastq_dir) exit 1, "basecall=false but --fastq_dir is not set"
    }
    if (!barcode_kit) exit 1, "Barcode kit not set. Provide --barcode_kit or --barcode_kit_name."

    Channel.fromPath(params.reference, checkIfExists: true).set { reference }
    Channel.fromPath(params.target_regions_bed, checkIfExists: true).set { target_regions_bed }

    Channel.fromPath(params.additional_metadata, checkIfExists: true)
        .ifEmpty { exit 1, "${params.additional_metadata} appears to be empty or missing" }
        .splitCsv(header:true, sep:',')
        .map { m -> ["${m.barcode_kit}_${m.barcode}", m] }
        .set { additional_metadata }

    def long_reads_ch
    def pycoqc_json_ch

    if (do_basecall) {
        def raw_reads = Channel.fromPath("${params.raw_read_dir}/*.{fast5,pod5}", checkIfExists: true)
        BASECALLING(raw_reads, additional_metadata)
        long_reads_ch  = BASECALLING.out.long_reads_ch
        pycoqc_json_ch = BASECALLING.out.pycoqc_json
    } else {
        long_reads_ch = Channel
            .fromPath("${params.fastq_dir}/*/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
            .map { f ->
                def barcode = f.parent.name
                def full_barcode = "${barcode_kit}_${barcode}"
                tuple(full_barcode, f)
            }
            .groupTuple()
            .map { full_barcode, files ->
                def barcode_raw = full_barcode.replaceFirst("^${barcode_kit}_", "")
                def meta = [
                    barcode_kit: barcode_kit,
                    barcode    : barcode_raw,
                    sample_id  : full_barcode,
                    ID         : full_barcode,
                    id         : full_barcode
                ]
                tuple(meta, files.sort())
            }
            .ifEmpty { exit 1, "No FASTQ files found under ${params.fastq_dir}" }
            | CAT_FASTQ

        pycoqc_json_ch = Channel.empty()
    }

    INDEX_REF(reference) | set { reference_index_ch }

    PRE_MAP_QC_PRE_TRIM(long_reads_ch)

    long_reads_ch
        .filter { meta, reads -> (meta.barcode ?: "").toString().toLowerCase() != "unclassified" }
        .set { remove_unclassified_for_mapping }

    PROCESS_FILTER_READS(remove_unclassified_for_mapping)
    PRE_MAP_QC_POST_TRIM(PROCESS_FILTER_READS.out.trimmed_reads)
    MAPPING(reference, PROCESS_FILTER_READS.out.trimmed_reads)

    FILTER_BAM(MAPPING.out.mapped_reads_bam, target_regions_bed)

    POST_FILTER_QC(
        FILTER_BAM.out.on_target_reads_bam,
        FILTER_BAM.out.off_target_reads_bam,
        target_regions_bed
    )

    def do_variant_calling = (params.clair3_model ?: "").toString().trim()

    if (do_variant_calling) {
        CALL_VARIANTS(
            FILTER_BAM.out.on_target_reads_bam.combine(target_regions_bed),
            reference_index_ch
        )
    } else {
        log.warn "Skipping CALL_VARIANTS: --clair3_model not provided"
    }

    MULTIQC(
        pycoqc_json_ch.ifEmpty([]),
        PRE_MAP_QC_PRE_TRIM.out.ch_fastqc_raw_zip.collect { it[1] }.ifEmpty([]),
        PRE_MAP_QC_POST_TRIM.out.ch_fastqc_raw_zip.collect { it[1] }.ifEmpty([]),
        MAPPING.out.ch_samtools_stats.collect { it[1,2] }.ifEmpty([]),
        POST_FILTER_QC.out.ch_samtools_stats.collect { it[1,2] }.ifEmpty([])
    )

    log.info "Pipeline complete"
}