process CLAIR3_CALL {
    label "cpu_8"
    label "mem_32"
    label "time_12"

    container 'hkubal/clair3:v1.0.9'

    // Avoid saveAs closures referencing input-scope vars in directive context
    publishDir "${params.outdir}/variants/", mode: 'copy', overwrite: true, pattern: 'merge_output.gvcf.gz'
    publishDir "${params.outdir}/variants/", mode: 'copy', overwrite: true, pattern: 'merge_output.vcf.gz'
    publishDir "${params.outdir}/variants/logs/", mode: 'copy', overwrite: true, pattern: 'run_clair3.log'

    input:
    tuple val(meta), path(filtered_bam), path(bam_index), path(target_regions_bed), path(reference), path(reference_index), path(clair3_model)

    output:
    tuple val(meta), path("merge_output.gvcf.gz"), emit: clair3_gvcf_out
    tuple val(meta), path("merge_output.gvcf.gz"), path("merge_output.gvcf.gz.tbi"), path(reference), path(reference_index), emit: clair3_gvcf_ref_idx_ch
    tuple val(meta), path("merge_output.vcf.gz")
    path("run_clair3.log")

    script:
    """
    MODEL_PATH=\$(readlink -f ${clair3_model})

    run_clair3.sh \\
      --bam_fn=${filtered_bam} \\
      --ref_fn=${reference} \\
      --threads=${task.cpus} \\
      --platform="ont" \\
      --model_path="\${MODEL_PATH}" \\
      --output=. \\
      --sample_name=${meta.ID} \\
      --bed_fn=${target_regions_bed} \\
      --include_all_ctgs \\
      --haploid_precise \\
      --min_coverage=${params.clair3_min_coverage} \\
      --call_snp_only \\
      --print_ref_calls \\
      --gvcf \\
      --var_pct_full=1 \\
      --ref_pct_full=1 \\
      --var_pct_phasing=1 \\
      --no_phasing_for_fa
    """
}