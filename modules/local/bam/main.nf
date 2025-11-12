process REMOVE_DUPLICATES_SAM {
    tag { sample }
    label "mediumcpu"

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
    label "highcpu"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/${sample}/aligned", pattern: "*.bam", mode: 'copy'

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
    label "highcpu"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(sam_list)

    output:
    tuple val(sample), path("merged.sorted.sam")

    script:
    def sam_files = sam_list instanceof List ? sam_list : [sam_list]
    def sam_inputs = sam_files.collect { it.toString() }.join(' ')
    """
    samtools merge -@ ${task.cpus} -n merged.sam ${sam_inputs} | \\
        samtools sort -@ ${task.cpus} -n -o merged.sorted.sam -
    """
}
