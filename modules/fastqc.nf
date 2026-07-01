process FASTQC {
    label 'cpu_2'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/qc/fastqc_pre_trim", mode: 'copy', overwrite: true, pattern: "*_fastqc.{html,zip}"
    publishDir "${params.outdir}/qc/fastqc_post_trim", mode: 'copy', overwrite: true, pattern: "*_trimmed_fastqc.{html,zip}"

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc \\
        -f fastq \\
        --threads ${task.cpus} \\
        ${reads}
    """
}