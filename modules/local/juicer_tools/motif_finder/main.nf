// Process: APA (aggregated peak analysis)
process MOTIF_FINDER {
    tag "${sample}"
    label "gpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

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
            ${params.genome_id} \\
            ${motif_dir} \\
            ${loops_txt}
    fi
    """
}
