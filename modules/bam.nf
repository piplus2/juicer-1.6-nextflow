process REMOVE_DUPLICATES_SAM {
    tag { sample }
    label "hpc"
    label "samtools"

    publishDir "${params.output_dir}/${sample}/aligned", pattern: "*.bam"

    input:
    tuple val(sample), path(merged_nodups), path(merged_sorted_sam)

    output:
    path "merged_nodups.bam", emit: alignable

    script:
    """
    filter_sam_by_readname.awk ${merged_nodups} ${merged_sorted_sam} | \\
        samtools view -@ ${task.cpus} -bS - > merged_nodups.bam
    """
}

process MERGE_SORT_SAM {
    tag { sample }
    label "hpc"
    label "samtools"

    input:
    tuple val(sample), path(sam_list)

    output:
    path "merged.sorted.sam", emit: merged_sorted_bam

    script:
    """
    samtools merge -@ ${task.cpus} -n merged.sam ${sam_list}
    samtools sort -@ ${task.cpus} -n -o merged.sorted.sam merged.sam
    """
}
