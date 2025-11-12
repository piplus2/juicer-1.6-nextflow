process BLACKLIST_CHIMERIC {
    tag { "${sample}-${name}" }
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(aligned_sam), path(norm_res_txt_initial)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}_norm.txt"), path("${name}${params.ext}_abnorm.sam"), path("${name}${params.ext}_unmapped.sam"), path("${name}${params.ext}_norm.sam"), path("${name}${params.ext}_norm.txt.res.txt")

    script:
    def output_norm = "${name}${params.ext}_norm.txt"
    def output_abnorm = "${name}${params.ext}_abnorm.sam"
    def output_unmapped = "${name}${params.ext}_unmapped.sam"
    def output_norm_sam = "${name}${params.ext}_norm.sam"
    """
    touch ${output_abnorm} ${output_unmapped}

    cp ${norm_res_txt_initial} ${name}${params.ext}_norm.txt.res.txt

    chimeric_blacklist.awk -v "fname1"=${output_norm} \\
        -v "fname2"=${output_abnorm} \\
        -v "fname3"=${output_unmapped} \\
        -v "fname_sam"=${output_norm_sam} \\
        ${aligned_sam}
    """
}
