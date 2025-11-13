process COUNT_LIGATIONS {
    tag "${sample}-${name}"
    label "highmemory"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    tuple val(sample), val(name), path(output_norm_txt), path(output_linecount_txt)

    script:
    output_norm_txt = "${name}${params.ext}_initial_norm.txt.res.txt"
    output_linecount_txt = "${name}${params.ext}_linecount.txt"
    """
    export LC_ALL=C
    export LC_COLLATE=C

    usegzip=0
    if [[ "${read1}" == *.gz ]]; then
        usegzip=1
    fi

    if [[ -z "${params.ligation}" || "${params.ligation}" == "XXXX" ]]; then
        echo "Skipping ligation match counting due to unset or 'XXXX' ligation motif"
        num1=0
    else
        if [[ \$usegzip -eq 1 ]]; then
            num1=\$(paste <(gunzip -c ${read1}) <(gunzip -c ${read2}) | awk '!((NR+2)%4)' | grep -cE ${params.ligation})
        else
            num1=\$(paste ${read1} ${read2} | awk '!((NR+2)%4)' | grep -cE "${params.ligation}")
        fi
    fi

    if [[ \$usegzip -eq 1 ]]; then
        num2=\$(gunzip -c ${read1} | wc -l | awk '{print \$1}')
    else
        num2=\$(wc -l ${read1} | awk '{print \$1}')
    fi
    echo -ne "\${num1} " > ${output_norm_txt}
    echo "\${num2}" > ${output_linecount_txt}
    """
}
