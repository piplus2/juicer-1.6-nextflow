process BLACKLIST_CHIMERIC {
    tag { "${sample}-${name}" }
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(aligned_sam), path(norm_res_txt_initial)

    output:
    tuple val(sample), val(name), path(output_norm), path(output_abnorm), path(output_unmapped), path(output_norm_sam), path(output_norm_res_txt)

    script:
    output_norm = "${name}${params.ext}_norm.txt"
    output_abnorm = "${name}${params.ext}_abnorm.sam"
    output_unmapped = "${name}${params.ext}_unmapped.sam"
    output_norm_sam = "${name}${params.ext}_norm.sam"
    output_norm_res_txt = "${name}${params.ext}_norm.txt.res.txt"
    """
    touch ${output_abnorm} ${output_unmapped}

    cp ${norm_res_txt_initial} ${output_norm_res_txt}
    chimeric_blacklist.awk -v "fname1"=${output_norm} \\
        -v "fname2"=${output_abnorm} \\
        -v "fname3"=${output_unmapped} \\
        -v "fname_sam"=${output_norm_sam} \\
        ${aligned_sam}
    """
}
