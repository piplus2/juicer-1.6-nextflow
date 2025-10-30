process REMOVE_DUPLICATES_SAM {
    tag { sample }
    label "hpc"

    input:
    tuple val(sample), path(merged_nodups), path(merged_sorted_sam)

    output:
    tuple val(sample), path("dedup.sam")

    script:
    """
    filter_sam_by_readname.awk ${merged_nodups} ${merged_sorted_sam} > dedup.sam
    """
}

process SAM_TO_BAM {
    tag { sample }
    label "hpc"
    label "samtools"

    publishDir "${params.output_dir}/${sample}/aligned", pattern: "*.bam", mode: 'copy'

    input:
    tuple val(sample), path(sam_file)

    output:
    tuple val(sample), path("${sam_file.baseName}.bam")

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam_file} -o ${sam_file.baseName}.bam
    """
}

process MERGE_SORT_SAM {
    tag { sample }
    label "hpc"
    label "samtools"

    input:
    tuple val(sample), path(sam_list)

    output:
    tuple val(sample), path("merged.sorted.sam")

    script:
    """
    samtools merge -@ ${task.cpus} -n merged.sam ${sam_list}
    samtools sort -@ ${task.cpus} -n -o merged.sorted.sam merged.sam
    """
}
