nextflow.enable.dsl = 2

include { process_fragments     } from './subworkflows/align.nf'
include { postprocessing        } from './subworkflows/postprocessing.nf'
include { REMOVE_DUPLICATES_SAM } from './modules/bam.nf'
include { GEN_HIC_FILES         } from './modules/hic.nf'
include { STATS                 } from './modules/stats.nf'
include { MAKE_HEADERFILE       } from './modules/header.nf'
include { MERGE_SORT ; REMOVE_DUPLICATES } from './modules/fragments.nf'
include { MERGE_SORT_SAM        } from './modules/bam.nf'

workflow {
    // Get the pairs of fastq files for all samples
    fastq_pairs = Channel.fromPath("${params.input_dir}/fastq/*_R{1,2}_001.fastq.gz")
        .map { file ->
            def sample = file.parent.parent.name
            def pair_id = file.name.replaceAll(/_R[12]_001\.fastq\.gz$/, "")
            tuple(sample, pair_id, file)
        }
        .groupTuple(by: [0, 1])
        .map { sample, pair_id, file ->
            def read1 = file.find { it.name.contains("_R1_") }
            def read2 = file.find { it.name.contains("_R2_") }
            tuple(sample, pair_id, read1, read2)
        }


    // Create one header file for each sample
    out_header = MAKE_HEADERFILE(fastq_pairs.map { it -> it[0] }.distinct())

    /* --------------- Process sample fragments (paired fastq) -------------- */

    frag_results = process_fragments(fastq_pairs)

    chimeric_reads = frag_results.chimeric_output
    sorted_fragments = frag_results.sorted_fragments

    // ( sample, sort_list ) -> ( sample, merged_sort_txt )
    sort_files_by_sample = sorted_fragments.groupTuple(by: 0)

    /* ------------------- Merge sorted and aligned files ------------------- */

    merged_sort = MERGE_SORT(sort_files_by_sample)

    /* -------------------------- Remove duplicates ------------------------- */

    merged_nodups = REMOVE_DUPLICATES(merged_sort)
    // Take tuple (sample, merged_nodups_txt) to process
    nodups = merged_nodups.map { it -> it[0..1] }

    /* ----------------------------- Process BAM ---------------------------- */

    // Merge the normal SAM files from the chimeric step for each sample
    norm_sam_by_sample = chimeric_reads
        .map { it -> tuple(it[0], it[2]) }
        .groupTuple(by: 0)
        .map { sample, sams -> tuple(sample, sams.join(" ")) }

    merged_sam = MERGE_SORT_SAM(norm_sam_by_sample)
    // Join (sample, merged_nodups_txt) with (sample, merged_sorted_sam)
    merged_nodups_by_sample = nodups.join(merged_sam)

    REMOVE_DUPLICATES_SAM(merged_nodups_by_sample)

    /* ----------------------------- Statistics ----------------------------- */

    // Remove name from chimeric tuple and group by sample id
    chimeric_by_sample = chimeric_reads
        .map { it -> tuple(it[0], it[3], it[4], it[5]) }
        .groupTuple(by: 0)

    stats_input = out_header
        .join(chimeric_by_sample)
        .join(nodups)

    stats_output = STATS(stats_input)

    /* -------------------------- Create HIC files -------------------------- */

    // (sample, inter, inter_30, inter_hists,  merged_nodups_txt) ->
    // (sample, inter.hic, inter_30.hic, inter_hists_30)

    hic_input = stats_output
        .map { it -> it[0..3] }
        .join(nodups, failOnMismatch: true)

    GEN_HIC_FILES(hic_input)

    inter_30_hic = GEN_HIC_FILES.out.map { it -> tuple(it[0], it[2]) }

    /* --------------------------- Postprocessing --------------------------- */

    postprocessing(inter_30_hic)
}
