process MAKE_HEADERFILE {
    tag "${sample}"
    label "smallcpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    val sample

    output:
    tuple val(sample), path(header_file)

    script:
    header_file = "header"
    """
    date > ${header_file}
    echo -ne 'Experiment description: ${params.about};' >> ${header_file}
    echo -ne "Juicer version: ${params.juicer_version};" >> ${header_file}
    bwa 2>&1 | awk '\$1=="Version:"{printf(" BWA %s; ", \$2)}' >> ${header_file}
    echo -ne "${params.threads} threads; " >> ${header_file}
    java -version 2>&1 | awk 'NR==1{printf("%s; ", \$0);}' >> ${header_file}
    juicer_tools -V 2>&1 | awk '\$1=="Juicer" && \$2=="Tools"{printf("%s; ", \$0);}' >> ${header_file}
    echo "" >> ${header_file}
    """
}
