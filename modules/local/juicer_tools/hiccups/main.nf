// Process: HiCCUPS (loop calling, requires GPU)
process HICCUPS {
    tag "${sample}"
    label "gpu"
    label "juicertools_1_22"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "${inter_30_hic.simpleName}_loops"

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    tuple val(sample), path(inter_30_hic), path(out_loops_dir, type: 'dir')

    script:
    out_loops_dir = "${inter_30_hic.simpleName}_loops"
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    mkdir -p ${out_loops_dir}

    if [[ ${params.use_gpu} == false ]]; then
        juicer_tools hiccups --cpu --threads ${task.cpus} \\
            ${inter_30_hic} \\
            ${out_loops_dir}
    else
        if hash nvcc 2> /dev/null; then
            juicer_tools hiccups --threads 0 \\
                ${inter_30_hic} \\
                ${out_loops_dir}
        else
            echo "ERROR: GPUs are required for HiCCUPS, but CUDA is not available." >&2
            exit 1
        fi
    fi
    """
}
