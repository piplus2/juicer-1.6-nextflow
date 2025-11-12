process BWA_ALIGN {
    tag "${sample}-${name}"
    label "highcpu"
    label "bwa"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}.sam")

    script:
    def output_aligned = "${name}${params.ext}.sam"
    """
    bwa-mem2 mem -SP5M -t ${task.cpus} ${params.reference} \\
        ${read1} ${read2} > ${output_aligned}
    """
}
