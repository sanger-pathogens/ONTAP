process GUNZIP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'ubuntu:20.04'

    input:
    tuple val(meta), path(zipped_gvcf), path(gvcf_index), path(reference), path(reference_index)

    output:
    tuple val(meta), path("*.gvcf"), path(gvcf_index), path(reference), path(reference_index), emit: gunzip

    script:
    """
    gunzip -c ${zipped_gvcf} > ${meta.ID}.gvcf
    """
}
