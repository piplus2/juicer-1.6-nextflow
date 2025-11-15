// Process: HiCCUPS (loop calling, requires GPU)
process HICCUPS {
    tag "${sample}"
    label "gpu"

    container 'docker://pinglese6022/juicer_tools:1.22.01'

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "${inter_30_hic.simpleName}_loops"

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    tuple val(sample), path(inter_30_hic), path(out_loops_dir, type: 'dir')

    script:
    out_loops_dir = "${inter_30_hic.simpleName}_loops"
    """
    export LC_ALL=en_US.UTF-8

    mkdir -p ${out_loops_dir}
    if hash nvcc 2> /dev/null; then
        /usr/local/bin/juicer_tools -Xmx${params.java_mem} hiccups --threads 0 --ignore-sparsity \\
            ${inter_30_hic} \\
            ${out_loops_dir}
    else
        echo "ERROR: GPUs are required for HiCCUPS, but CUDA is not available." >&2
        exit 1
    fi
    """
}
