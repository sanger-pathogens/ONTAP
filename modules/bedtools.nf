process BEDTOOLS_GENOMECOV {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/qc/coverage/bedtools_genome_coverage", mode: 'copy', overwrite: true, pattern: "*.bedGraph"

    conda "bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bam_index)
    val(qc_stage)

    output:
    path("*.bedGraph"), emit: genome_cov_ch

    script:
    """
    bedtools genomecov -split -ibam ${sorted_bam} -bga | \
    bedtools sort > ${meta.ID}.bedGraph
    """
}

process BEDTOOLS_COVERAGE {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/qc/coverage/bedtools_coverage", mode: 'copy', overwrite: true, pattern: "*coverage.bed"

    conda "bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bam_index), path(target_regions_bed)
    val(qc_stage)

    output:
    path("*.bed"), emit: genome_cov_ch

    script:
    """
    bedtools coverage -hist -abam ${sorted_bam} -b ${target_regions_bed} > ${meta.ID}_coverage.bed
    """
}