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
