process REMOVE_DUPLICATES_SAM {
    tag { sample }
    label "mediumcpu"

    input:
    tuple val(sample), path(merged_nodups), path(merged_sorted_sam)

    output:
    tuple val(sample), path(output_sam)

    script:
    output_sam = "dedup.sam"
    """
    filter_sam_by_readname.awk ${merged_nodups} ${merged_sorted_sam} > ${output_sam}
    """
}
