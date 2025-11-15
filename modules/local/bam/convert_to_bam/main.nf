process SAM_TO_BAM {
    tag { sample }
    label "highcpu"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'biocontainers/samtools:1.21--h50ea8bc_0'}"

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
