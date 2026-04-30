include { SAMTOOLS_SORT         } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_MERGE        } from '../../../modules/nf-core/samtools/merge/main.nf'
include { PICARD_FILTERSAMREADS } from '../../../modules/nf-core/picard/filtersamreads/main.nf'

workflow EXPORT_BAM {
    take:
    aligned_sams // Expect a tuple of (sample, name, norm_txt, abnorm_sam, unmapped_sam, norm_sam, norm_res_txt)
    dedup        // Expect a tuple of (sample, merged_nodups, dups, opt_dups)

    main:
    // Group the normalized SAM files by sample and create the nf-core [meta, [files]] structure
    // We use groupTuple to group by sample, and then map to create the tuple structure expected by SAMTOOLS_VIEW
    ch_to_sort = aligned_sams
        .map { sample, _name, sam ->
            def meta = [id: sample, sample: sample]
            return [meta, sam]
        }
        .groupTuple(by: 0)

    // Sort (and implicitly merge)
    // SAMTOOLS_SORT can take a list of SAM files.
    // It will output a single merged, sorted BAM file.
    SAMTOOLS_SORT(ch_to_sort, [[:], [], []], 'bai')

    // Get the deduplicated IDs from merged_nodups
    ch_id_list = dedup.map { sample, nodups ->
        def id_list = file("${sample}_allowed_ids.txt")
        // Column 15 and 16 contain the QNAME in merged_nodups
        def text = ""
        nodups.splitEachLine("\\s+") { cols ->
            if (cols.size() >= 16) {
                text += cols[14] + "\n" + cols[15] + "\n"
            }
        }
        id_list.text = text
        return [[id: sample, sample: sample], id_list]
    }

    // Filter the merged BAM to keep only reads with QNAME in the allowed IDs list
    PICARD_FILTERSAMREADS(SAMTOOLS_SORT.out.bam, ch_id_list.map { _meta, ids -> ids })

    emit:
    bam = PICARD_FILTERSAMREADS.out.bam
}
