// Process: Arrowhead (contact domain calling)
process ARROWHEAD {
    tag "${sample}"
    label "gpu"

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    path "${inter_30_hic.simpleName}_contact_domains", type: 'dir'

    script:
    def contact_domains = "${inter_30_hic.simpleName}_contact_domains"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8
    export JAVA_HOME=/work/pinglese/tools/jdk-17.0.14
    export PATH=\${JAVA_HOME}/bin:\${PATH}

    mkdir -p ${contact_domains}

    juicer_tools arrowhead \\
        --threads 0 \\
        --ignore-sparsity \\
        ${inter_30_hic} \\
        ${contact_domains}
    """
}


// Process: HiCCUPS (loop calling, requires GPU)
process HICCUPS {
    tag "${sample}"
    label "gpu"

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy', pattern: "${inter_30_hic.simpleName}_loops"

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    tuple val(sample), path(inter_30_hic), path("${inter_30_hic.simpleName}_loops", type: 'dir')

    script:
    def out_loops = "${inter_30_hic.simpleName}_loops"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8
    export JAVA_HOME=/work/pinglese/tools/jdk-17.0.14
    export PATH=\${JAVA_HOME}/bin:\${PATH}

    mkdir -p ${out_loops}

    if hash nvcc 2> /dev/null; then
        juicer_tools hiccups \\
        --threads 0 \\
        --ignore-sparsity \\
        ${resolutions} \\
        ${inter_30_hic} \\
        ${out_loops} \\
    else
        echo "ERROR: GPUs are required for HiCCUPS, but CUDA is not available." >&2
        exit 1
    fi
    """
}


// Process: APA (aggregated peak analysis)
process MOTIF_FINDER {
    tag "${sample}"
    label "gpu"

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic), path(merged_loops_dir)

    output:
    path "apa_results", type: 'dir'
    path "${inter_30_hic.simpleName}_loops_with_motifs.bedpe"

    script:
    def motif_dir = params.motif_dir == null || params.motif_dir.toString() == "" ? "" : "${params.motif_dir}"
    def loops_txt = "${inter_30_hic.simpleName}_loops.txt"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    export JAVA_HOME=/work/pinglese/tools/jdk-17.0.14
    export PATH=\${JAVA_HOME}/bin:\${PATH}
    mkdir -p "apa_results"

    cp ${merged_loops_dir}/merged_loops.bedpe ${loops_txt}

    juicer_tools apa \\
        --threads 1 \\
        ${inter_30_hic} \\
        ${merged_loops_dir}/merged_loops.bedpe \\
        "apa_results"

    if [[ -z "${motif_dir}" ]]; then
        echo "***! Can't find folder ${motif_dir}"
        echo "***! Not running motif finder"
        touch ${loops_txt}
    else
        juicer_tools motifs \\
            ${params.genomeID} \\
            ${motif_dir} \\
            ${loops_txt} \\
    fi
    """
}

workflow postprocessing {
    take:
    inter_30

    main:

    ARROWHEAD(inter_30)

    loops_dir = HICCUPS(inter_30)
    // Keep the loops dirs that contain the file merged_loops.bedpe
    good_loops_dir = loops_dir.filter { _sample, _inter_30, _loops_dir ->
        file("${_loops_dir}/merged_loops.bedpe").exists()
    }

    MOTIF_FINDER(good_loops_dir)

    emit:
    null
}
