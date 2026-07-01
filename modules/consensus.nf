process CURATE_CONSENSUS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/curated_consensus", pattern: "*.fasta", mode: 'copy', overwrite: true

    conda 'conda-forge::python=3.10.2'
    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), file(gvcf_final), file(gvcf_index), path(reference), path(ref_index), path(target_regions_bed)

    output:
    tuple val(meta), path("${meta.ID}.fasta"),  emit: multi_fasta
    tuple val(meta), path("${meta.ID}_multi_locus.fasta"), emit: full_consensus
    tuple val(meta), path("${meta.ID}_wg.fasta"), emit: wg_consensus
    //path("*.fa"), emit: per_loci

    script:
    """
    gvcf_to_fasta.py -r ${reference} -v ${gvcf_final} --fasta_id ${meta.ID} -b ${target_regions_bed} -n --multi_locus --multifasta --whole_genome_fasta -o ${meta.ID} --min_ref_gt_qual ${params.min_ref_gt_qual} --min_alt_gt_qual ${params.min_alt_gt_qual}
    # below outputs per-sample, per-locus fasta files THAT HAVE REF BASE AS DEFAULT BASE ON POSITIONS THAT ARE NOT SUPPORTED INSTEAD OF AN N CHAR
    # gvcf_to_fasta.py -r ${reference} -v ${gvcf_final} --fasta_id ${meta.ID} -b ${target_regions_bed} --singlefasta -o ${meta.ID} --min_ref_gt_qual ${params.min_ref_gt_qual} --min_alt_gt_qual ${params.min_alt_gt_qual}
    """
}
