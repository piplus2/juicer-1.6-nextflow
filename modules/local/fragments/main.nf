process CONVERT_FRAGMENTS {
    tag "${sample}-${name}"
    label "smallcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(norm_txt)

    output:
    tuple val(sample), val(name), path(output_frag_txt)

    script:
    output_frag_txt = "${name}${params.ext}.frag.txt"
    """
    if [ "${params.site}" != "none" ] && [ -e "${params.site_file}" ]; then
        fragment.pl \\
            ${norm_txt} \\
            ${output_frag_txt} \\
            ${params.site_file}
    elif [ "${params.site}" == "none" ] || [ "${params.nofrag}" -eq 1 ]; then
        awk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\$i);}printf("\n");}' ${norm_txt} > ${output_frag_txt}
    else
        echo "***! No ${norm_txt} file created"
        exit 1
    fi
    """
}
