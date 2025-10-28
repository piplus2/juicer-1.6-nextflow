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
    // Discover paired FASTQ files once and reuse the resulting tuples downstream
    fastq_pairs = Channel.fromFilePairs(
            "${params.input_dir}/fastq/*_R{1,2}_001.fastq.gz",
            checkIfExists: true
        )
        .map { pair_id, reads ->
            def sorted_reads = reads.sort { it.name }
            assert sorted_reads.size() == 2 : "Expected two FASTQ mates for ${pair_id}"

            def sample = sorted_reads[0].parent.parent.name
            tuple(sample, pair_id, sorted_reads[0], sorted_reads[1])
        }

    // Split channel consumption between header generation and fragment processing
    def fastq_for_header = null
    def fastq_for_processing = null
    fastq_pairs.into {
        fastq_for_header
        fastq_for_processing
    }

    // Create one header file for each sample
    out_header = MAKE_HEADERFILE(
        fastq_for_header.map { sample, _pair_id, _read1, _read2 -> sample }.distinct()
    )

    /* --------------- Process sample fragments (paired fastq) -------------- */

    frag_results = process_fragments(fastq_for_processing)

    chimeric_reads = frag_results.chimeric_output
    sorted_fragments = frag_results.sorted_fragments

    // ( sample, sort_list ) -> ( sample, merged_sort_txt )
    sort_files_by_sample = sorted_fragments.groupTuple(by: 0)

    /* ------------------- Merge sorted and aligned files ------------------- */

    merged_sort = MERGE_SORT(sort_files_by_sample)

    /* -------------------------- Remove duplicates ------------------------- */

    merged_nodups = REMOVE_DUPLICATES(merged_sort)
    // Take tuple (sample, merged_nodups_txt) to process
    nodups = merged_nodups.map { sample, merged_nodups_txt, _dups_txt, _opt_dups_txt ->
        tuple(sample, merged_nodups_txt)
    }

    /* ----------------------------- Process BAM ---------------------------- */

    // Merge the normal SAM files from the chimeric step for each sample
    norm_sam_by_sample = chimeric_reads
        .map { sample, _name, _norm_txt, _abnorm_sam, _unmapped_sam, _norm_res_txt, norm_sam ->
            tuple(sample, norm_sam)
        }
        .groupTuple(by: 0)
        .map { sample, sams -> tuple(sample, sams.join(" ")) }

    merged_sam = MERGE_SORT_SAM(norm_sam_by_sample)
    // Join (sample, merged_nodups_txt) with (sample, merged_sorted_sam)
    merged_nodups_by_sample = nodups.join(merged_sam)

    REMOVE_DUPLICATES_SAM(merged_nodups_by_sample)

    /* ----------------------------- Statistics ----------------------------- */

    // Remove name from chimeric tuple and group by sample id
    chimeric_by_sample = chimeric_reads
        .map { sample, _name, _norm_txt, abnorm_sam, unmapped_sam, norm_res_txt, _norm_sam ->
            tuple(sample, norm_res_txt, abnorm_sam, unmapped_sam)
        }
        .groupTuple(by: 0)

    // out_header: (sample, header_file)
    // chimeric_by_sample: (sample, [norm_res_txt, abnorm_sam, unmapped_sam] )
    // nodups: (sample, merged_nodups_txt)
    stats_input = out_header
        .join(chimeric_by_sample)
        .join(nodups)

    // stats_output: (sample, inter, inter_30, inter_hists, collisions, abnorm_sam, unmapped_sam)
    stats_output = STATS(stats_input)

    /* -------------------------- Create HIC files -------------------------- */

    // (sample, inter, inter_30, inter_hists,  merged_nodups_txt) ->
    // (sample, inter.hic, inter_30.hic, inter_hists_30)

    hic_input = stats_output
        .map { sample, inter, inter_30, inter_hists, _collisions, _abnorm_sam, _unmapped_sam ->
            tuple(sample, inter, inter_30, inter_hists)
        }
        .join(nodups, failOnMismatch: true)

    GEN_HIC_FILES(hic_input)

    inter_30_hic = GEN_HIC_FILES.out.map { sample, _inter_hic, inter_30_hic_file, _inter_30_hists ->
        tuple(sample, inter_30_hic_file)
    }

    /* --------------------------- Postprocessing --------------------------- */

    postprocessing(inter_30_hic)
}
