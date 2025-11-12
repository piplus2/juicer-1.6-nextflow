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
