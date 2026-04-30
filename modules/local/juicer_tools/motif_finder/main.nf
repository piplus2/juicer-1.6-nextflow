// Process: APA (aggregated peak analysis)
process MOTIF_FINDER {
    tag "${sample}"
    label "highcpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic), path(merged_loops_dir)
    path motif_dir

    output:
    path "apa_results", type: 'dir'
    path "${inter_30_hic.simpleName}_loops_with_motifs.bedpe"

    script:
    def loops_txt = "${inter_30_hic.simpleName}_loops.txt"
    def motif_dir_val = motif_dir ? motif_dir : ""
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    mkdir -p "apa_results"

    cp ${merged_loops_dir}/merged_loops.bedpe ${loops_txt}

    juicer_tools apa \\
        --threads 1 \\
        ${inter_30_hic} \\
        ${merged_loops_dir}/merged_loops.bedpe \\
        "apa_results"

    if [[ -z "${motif_dir_val}" ]]; then
        log.error("Motif directory not provided. Skipping motif finding.")
        touch ${loops_txt}
    else if [[ ! -d "${motif_dir_val}" ]]; then
        log.error("Motif directory '${motif_dir_val}' does not exist. Skipping motif finding.")
        touch ${loops_txt}
    else
        juicer_tools motifs \\
            ${params.genome_id} \\
            ${motif_dir_val} \\
            ${loops_txt}
    fi
    """
}
