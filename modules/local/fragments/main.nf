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
    // Correctly convert the boolean parameter to a string for Bash
    def nofrag_val = params.nofrag ? "true" : "false"
    """
    if [ "${params.site}" != "none" ] && [ -f "${params.site_file}" ]; then
        echo "Restriction enzyme site provided, assigning fragments based on site cutting"
        fragment.pl \\
            ${norm_txt} \\
            ${output_frag_txt} \\
            ${params.site_file}

    elif [ "${params.site}" == "none" ] || [ "${nofrag_val}" == "true" ]; then
        echo "No restriction enzyme site provided or nofrag option set to true, skipping fragment assignment"
        # Escaping \$ for Gawk and using \\n for the newline character
        gawk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\$i);}printf("\\n");}' ${norm_txt} > ${output_frag_txt}

    else
        echo "Error: Configuration mismatch. Site is ${params.site} but site file not found, or internal logic failed."
        exit 1
    fi
    """
}
