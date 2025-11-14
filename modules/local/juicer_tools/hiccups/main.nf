// Process: HiCCUPS (loop calling, requires GPU)
process HICCUPS {
    tag "${sample}"
    label "gpu"

    container 'docker://pinglese6022/juicer_tools:1.22.01'

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "${inter_30_hic.simpleName}_loops"

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    tuple val(sample), path(inter_30_hic), path("${inter_30_hic.simpleName}_loops", type: 'dir')

    script:
    def out_loops = "${inter_30_hic.simpleName}_loops"
    """
    export LC_ALL=en_US.UTF-8

    mkdir -p ${out_loops}

    if hash nvcc 2> /dev/null; then
        args="--threads 0 --ignore-sparsity"
        if [[ -n "${params.resolutions}" ]]; then
            args+="-r \"${params.resolutions}\""
        fi
        args+="\"${inter_30_hic}\" \"${out_loops}\""

        /usr/local/bin/juicer_tools -Xmx${params.java_mem} hiccups ${args}
    else
        echo "ERROR: GPUs are required for HiCCUPS, but CUDA is not available." >&2
        exit 1
    fi
    """
}
