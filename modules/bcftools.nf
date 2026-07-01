process MERGE_GVCF {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/variants/merged_gvcf/", mode: 'copy', overwrite: true, pattern: "*_merged.vcf.gz"

    conda "bioconda::bcftools=1.20"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    path("variant_vcfs*.gvcf.gz")
    path(target_regions_bed)

    output:
    path(out_vcf), emit: merged_vcf

    script:
    //workflow start is ugly 2024-02-29T12:01:26.233465Z split on T for time and take just the date
    date = "${workflow.start}".split('T')[0]

    out_vcf = "${workflow.runName}_${date}_merged.vcf.gz"
    """
    find . -name "*.gvcf.gz" -exec bcftools index "{}" \\;

    bcftools merge -R ${target_regions_bed} -0 -o ${out_vcf} -O z -f PASS *.gvcf.gz
    """
}

process BCFTOOLS_QUERY {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/variants/merged_gvcf/", mode: 'copy', overwrite: true, pattern: "*.tsv"

    conda "bioconda::bcftools=1.20"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    path(merged_vcf)

    output:
    path(out_tsv), emit: parsed_tsv

    script:
    //workflow start is ugly 2024-02-29T12:01:26.233465Z split on T for time and take just the date
    date = "${workflow.start}".split('T')[0]

    out_tsv = "${workflow.runName}_${date}_merged.tsv"
    """
    bcftools query -f'[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%TGT\\t%GQ\\t%DP\\t%AD\\t%AF\\n]' ${merged_vcf} > ${out_tsv}
    """
}
