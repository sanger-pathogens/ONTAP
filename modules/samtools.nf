process CONVERT_TO_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir path: "${params.outdir}/fastqs/", enabled: (params.save_fastqs in [true, "true"]), mode: 'copy', overwrite: true, pattern: "*.fastq.gz"

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(reads_bam)

    output:
    tuple val(meta), path(fastq_output), emit: reads_fastq

    script:
    fastq_output = "${meta.ID}.fastq.gz"
    """
    samtools fastq -@ ${task.cpus} -0 ${fastq_output} ${reads_bam}
    """
}

process MERGE_BAMS_FOR_SUMMARY {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    path("*.bam")

    output:
    path(combined_bam), emit: summary_bam

    script:
    combined_bam = "merged.bam"
    """
    samtools merge -@ ${task.cpus} -o ${combined_bam} *.bam
    """
}

process MANAGE_DUPLICATES_FROM_BAMS {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam), path(duplicates_list)
    val(mode)

    output:
    tuple val(meta), path(final_bam), emit: bam

    script:
    final_bam = "${bam.simpleName}_clean.bam"
    if (mode == "keep")
        """
        samtools view -N ${duplicates_list} -o ${final_bam} ${bam}
        """
    else if (mode == "remove")
        """
        samtools view -N ^${duplicates_list} -o ${final_bam} ${bam}
        """
}

process CONVERT_TO_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(mapped_reads)

    output:
    tuple val(meta), path("${mapped_reads_bam}"), emit: mapped_reads_bam

    script:
    mapped_reads_bam = "${meta.ID}.bam"
    """
    samtools view \\
        -@ ${task.cpus} \\
        -bS -F 2308 \\
        --min-MQ 50 \\
        -o ${mapped_reads_bam} \\
        ${mapped_reads}
    """
}

process SAMTOOLS_SORT {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/mapped_reads", enabled: (params.keep_sorted_bam in [true, "true"]), mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(mapped_reads_bam)

    output:
    tuple val(meta), path("${sorted_reads}"), emit: sorted_reads

    script:
    sorted_reads = "${meta.ID}_sorted.bam"
    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${sorted_reads} \
        ${mapped_reads_bam}
    """
}

process REMOVE_OFF_TARGET_READS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(sorted_reads_bam), path(sorted_reads_bai), path(target_regions_bed)

    output:
    tuple val(meta), path("${on_target_reads_bam}"), emit: on_target_reads_bam
    tuple val(meta), path("${off_target_reads_bam}"), emit: off_target_reads_bam

    script:
    on_target_reads_bam = "${meta.ID}_on_target.bam"
    off_target_reads_bam = "${meta.ID}_off_target.bam"
    """
    samtools view \
        -@ ${task.cpus} \
        -L ${target_regions_bed} \
        -o ${on_target_reads_bam} \
        -U ${off_target_reads_bam} \
        ${sorted_reads_bam}
    """
}

process ON_AND_OFF_TARGET_STATS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(on_target_reads_bam), path(on_target_reads_bai), path(off_target_reads_bam), path(off_target_reads_bai)

    output:
    tuple val(meta), path("${on_and_off_target_stats}"), emit: on_and_off_target_stats

    script:
    on_and_off_target_stats = "${meta.ID}_on_and_off_target_stats.csv"
    """
    on_and_off_target_stats.sh ${meta.ID} ${on_target_reads_bam} ${off_target_reads_bam} ${on_and_off_target_stats}
    """
}

process INDEX_REF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/sorted_ref", mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    path(reference)

    output:
    tuple path(reference), path("${faidx}"), emit: ref_index

    script:
    faidx = "${reference}.fai"
    """
    samtools faidx "${reference}" > "${faidx}"
    """
}

process SAMTOOLS_INDEX_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    publishDir "${params.outdir}/mapped_reads", enabled: (params.keep_bam_files in [true, "true"]), mode: 'copy', overwrite: true, pattern: "*.bai"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path(bam_file), path("*.bai"), emit: bam_index

    script:
    """
    samtools index -@ ${task.cpus} *.bam
    """
}

process GET_READLENGTH_DISTRIBUTION {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/qc/readlengths", mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(sorted_bam), path(bam_index)
    val(qc_stage)

    output:
    tuple val(meta), path("*.read-lengths.tsv"), emit: readlengths

    script:
    """
    samtools view -@ ${task.cpus} ${sorted_bam} | \
    awk '{print length(\$10)}' | sort | uniq -c | sort -n -k 2 | awk -v OFS='\t' '{print \$2,\$1}' \
    > ${meta.ID}.read-lengths.tsv
    """
}

process SAMTOOLS_DEPTH {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/qc/coverage/samtools_depth", mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam_file), path(bam_index)
    val(qc_stage)

    output:
    tuple val(meta), path(coverage_report), emit: samtools_coverage

    script:
    coverage_report = "${meta.ID}_samtools_depth.tsv"
    """
    samtools depth -@ ${task.cpus} -aa *.bam -o ${coverage_report}
    """
}

process SAMTOOLS_STATS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/qc/samtools_stats", mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(mapped_reads_bam)
    val(qc_stage)

    output:
    tuple val(meta), path(stats_file), path(flagstats_file), emit: stats_ch

    script:
    stats_file = "${meta.ID}.stats"
    flagstats_file = "${meta.ID}.flagstats"
    """
    samtools stats -@ ${task.cpus} "${mapped_reads_bam}" > "${stats_file}"
    samtools flagstats -@ ${task.cpus} "${mapped_reads_bam}" > "${flagstats_file}"
    """
}