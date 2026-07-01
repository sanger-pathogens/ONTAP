include { CONVERT_FAST5_TO_POD5; MERGE_POD5 } from '../modules/pod5.nf'
include { MODEL_DOWNLOAD; BASECALL; DEMUX; DORADO_SUMMARY; UNASSIGNED_SUMMARY } from '../modules/dorado.nf'
include {
    CONVERT_TO_FASTQ;
    MERGE_BAMS_FOR_SUMMARY;
    MANAGE_DUPLICATES_FROM_BAMS as REMOVE_DUPLICATES_FROM_BAMS;
    MANAGE_DUPLICATES_FROM_BAMS as KEEP_DUPLICATES_FROM_BAMS
} from '../modules/samtools.nf'
include { SUMMARY_DUPLICATES } from '../modules/summary_duplicates.nf'
include { PYCOQC } from '../modules/pycoqc.nf'

def validateSingleFormat(listOfFormats) {
    if (listOfFormats.size() != 1) {
        error("Multiple signal filetypes ${listOfFormats} found in '${params.raw_read_dir}'. Please separate filetypes into distinct directories and process independently.")
    }
}

workflow BASECALLING {
    take:
    raw_read_signal_files
    additional_metadata

    main:
    raw_read_signal_files
        .map { raw_read_signal_file -> tuple(raw_read_signal_file.extension, raw_read_signal_file) }
        | groupTuple
        | set { input_formats }

    input_formats
        .map { format, files -> format }
        | collect
        | map { format -> validateSingleFormat(format) }

    input_formats
        | branch { format, files ->
            fast5: format == "fast5"; return files
            pod5 : format == "pod5";  return files
        }
        | set { raw_files }

    CONVERT_FAST5_TO_POD5(raw_files.fast5)

    MERGE_POD5(raw_files.pod5)
        | mix(CONVERT_FAST5_TO_POD5.out.pod5_ch)
        | set { pod5_ch }

    if ((params.basecall_model_path ?: "") != "" && file(params.basecall_model_path).exists()) {
        pod5_ch
            | combine(Channel.fromPath(params.basecall_model_path))
            | BASECALL
    } else {
        pod5_ch
            | MODEL_DOWNLOAD
            | BASECALL
    }

    DEMUX(BASECALL.out.called_channel, params.barcode_kit_name)
        | transpose
        | map { barcode_kit, long_read_bam ->
            def meta = [:]
            meta.barcode_kit = barcode_kit
            meta.barcode = long_read_bam.simpleName.contains("barcode")
                ? long_read_bam.simpleName.split("barcode")[1]
                : long_read_bam.simpleName
            tuple(meta, long_read_bam)
        }
        | set { barcode_bam_ch }

    LONG_READ_QC(barcode_bam_ch)

    LONG_READ_QC.out.pycoqc_json
        | set { pycoqc_json }

    if (params.barcode_kit_name.size() == 1) {
        barcode_bam_ch.set { bam_ch }
    } else {
        SUMMARY_DUPLICATES(LONG_READ_QC.out.summary_channel, "remove")
            | set { duplicate_classified_list }

        barcode_bam_ch
            .filter { meta, long_bam -> long_bam.name != "unclassified.bam" }
            | combine(duplicate_classified_list)
            | set { marked_for_removal_bam }

        REMOVE_DUPLICATES_FROM_BAMS(marked_for_removal_bam, "remove")

        SORT_UNCLASSIFIED(LONG_READ_QC.out.unclassified_ch)

        REMOVE_DUPLICATES_FROM_BAMS.out.bam
            | set { bam_ch }
    }

    bam_ch
        | map { meta, reads -> ["${meta.barcode_kit}_${meta.barcode}", meta, reads] }
        | set { bam_by_barcode_ch }

    additional_metadata
        | join(bam_by_barcode_ch)
        | map { barcodekit_barcode, meta1, meta2, reads -> [meta1 + meta2, reads] }
        | set { bam_with_metadata_ch }

    if (params.read_format == "fastq") {
        if (params.barcode_kit_name.size() == 1) {
            CONVERT_TO_FASTQ(bam_with_metadata_ch)
                | set { long_reads_ch }
        } else {
            bam_with_metadata_ch
                .mix(SORT_UNCLASSIFIED.out.cleaned_unclassified)
                | CONVERT_TO_FASTQ
                | set { long_reads_ch }
        }
    } else {
        if (params.barcode_kit_name.size() == 1) {
            bam_with_metadata_ch
                | set { long_reads_ch }
        } else {
            bam_with_metadata_ch
                .mix(SORT_UNCLASSIFIED.out.cleaned_unclassified)
                | set { long_reads_ch }
        }
    }

    emit:
    long_reads_ch   = long_reads_ch
    summary_channel = LONG_READ_QC.out.summary_channel
    pycoqc_json     = pycoqc_json
}

workflow LONG_READ_QC {
    take:
    barcode_bam_ch

    main:
    barcode_bam_ch
        .branch { metaMap, long_bam ->
            unclassified: metaMap.barcode == "unclassified"; return long_bam
            classified  : true;                              return long_bam
        }
        .set { summarise_channel }

    summarise_channel.classified
        .collect()
        | MERGE_BAMS_FOR_SUMMARY
        | DORADO_SUMMARY
        | set { summary_channel }

    PYCOQC(summary_channel)
    PYCOQC.out.json
        | set { pycoqc_json }

    if (params.barcode_kit_name.size() == 1) {
        summarise_channel.unclassified
            | set { unclassified_ch }
    } else {
        summarise_channel.unclassified
            .collect()
            | set { unclassified_ch }
    }

    emit:
    summary_channel = summary_channel
    unclassified_ch = unclassified_ch
    pycoqc_json     = pycoqc_json
}

workflow SORT_UNCLASSIFIED {
    take:
    unclassified

    main:
    MERGE_BAMS_FOR_SUMMARY(unclassified)
        | UNASSIGNED_SUMMARY
        | set { unclassified_summary }

    SUMMARY_DUPLICATES(unclassified_summary, "keep")
        | set { unclassified_duplicates }

    MERGE_BAMS_FOR_SUMMARY.out.summary_bam
        | combine(unclassified_duplicates)
        | map { long_bam, sequencing_summary ->
            def unclassified_meta = [:]
            unclassified_meta.barcode_kit = "Multiple"
            unclassified_meta.barcode = "Unclassified"
            unclassified_meta.ID = "Unclassified_reads"
            tuple(unclassified_meta, long_bam, sequencing_summary)
        }
        | set { unclassified_marked }

    KEEP_DUPLICATES_FROM_BAMS(unclassified_marked, "keep")
        | set { cleaned_unclassified }

    emit:
    cleaned_unclassified = cleaned_unclassified
}