process MERGE_SORT_SAM {
    tag { sample }
    label "highcpu"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(sample), path(sam_list)

    output:
    tuple val(sample), path("merged.sorted.sam")

    script:
    def sam_files = sam_list instanceof List ? sam_list : [sam_list]
    def sam_inputs = sam_files.collect { it.toString() }
    if (sam_files.size() == 1) {
        """
        samtools sort -@ ${task.cpus} -n -o merged.sorted.sam ${sam_inputs[0]}
        """
    } else {
        def sam_inputs_str = sam_inputs.join(' ')
        """
        samtools merge -@ ${task.cpus} -n merged.sam ${sam_inputs_str} | \\
            samtools sort -@ ${task.cpus} -n -o merged.sorted.sam -
        """
    }
}
