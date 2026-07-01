include { CLAIR3_CALL } from '../modules/clair3.nf'
include { GUNZIP } from '../modules/helper_process.nf'
include { CURATE_CONSENSUS } from '../modules/consensus.nf'
include { MERGE_GVCF; BCFTOOLS_QUERY } from '../modules/bcftools.nf'
include { CONSTRUCT_PHYLO } from '../assorted-sub-workflows/tree_build/tree_build.nf'

workflow CALL_VARIANTS {
    take:
    filtered_reads_bam_with_bed
    reference_index_ch

    main:
    def clair3_model_path = (params.clair3_model ?: "").toString().trim()
    def target_bed_path   = (params.target_regions_bed ?: "").toString().trim()

    if (!clair3_model_path) {
        error "CALL_VARIANTS: --clair3_model is required but was not provided"
    }
    if (!target_bed_path) {
        error "CALL_VARIANTS: --target_regions_bed is required but was not provided"
    }

    def clair3_model_ch = Channel.fromPath(clair3_model_path, checkIfExists: true)
    def target_bed_ch   = Channel.fromPath(target_bed_path, checkIfExists: true)

    filtered_reads_bam_with_bed
        | combine(reference_index_ch)
        | combine(clair3_model_ch)
        | CLAIR3_CALL

    GUNZIP(CLAIR3_CALL.out.clair3_gvcf_ref_idx_ch)
        | combine(target_bed_ch)
        | CURATE_CONSENSUS

    MERGE_GVCF(
        CLAIR3_CALL.out.clair3_gvcf_out.map { metadata, path -> path }.collect(),
        target_bed_ch
    )

    BCFTOOLS_QUERY(MERGE_GVCF.out.merged_vcf)

    CURATE_CONSENSUS.out.full_consensus.collectFile { meta, file -> ["merged.fasta", file] }
        | CONSTRUCT_PHYLO

    emit:
    clair3_gvcf_out = CLAIR3_CALL.out.clair3_gvcf_out
}