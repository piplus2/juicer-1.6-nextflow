// Process: APA (aggregated peak analysis)
process MOTIF_FINDER {
    tag "${sample}"
    label "gpu"

    container 'docker://pinglese6022/juicer_tools:1.13.01'

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic), path(merged_loops_dir)
    path motif_dir

    output:
    path "apa_results", type: 'dir'
    path "${inter_30_hic.simpleName}_loops_with_motifs.bedpe"

    script:
    def loops_txt = "${inter_30_hic.simpleName}_loops.txt"
    """
    export LC_ALL=en_US.UTF-8
    mkdir -p "apa_results"

    cp ${merged_loops_dir}/merged_loops.bedpe ${loops_txt}

    /usr/local/bin/juicer_tools -Xmx${params.java_mem} apa \\
        --threads 1 \\
        ${inter_30_hic} \\
        ${merged_loops_dir}/merged_loops.bedpe \\
        "apa_results"

    if [[ -z "${motif_dir}" ]]; then
        echo "***! Can't find folder ${motif_dir}"
        echo "***! Not running motif finder"
        touch ${loops_txt}
    else
        /usr/local/bin/juicer_tools -Xmx${params.java_mem} motifs \\
            ${params.genome_id} \\
            ${motif_dir} \\
            ${loops_txt}
    fi
    """
}
