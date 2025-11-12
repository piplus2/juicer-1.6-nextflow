nextflow.enable.dsl = 2

include { process_fragments     } from '../subworkflows/local/align.nf'
include { postprocessing        } from '../subworkflows/local/postprocessing.nf'
include { REMOVE_DUPLICATES_SAM ; SAM_TO_BAM } from '../modules/local/bam.nf'
include { GEN_HIC_FILES         } from '../modules/local/hic.nf'
include { STATS                 } from '../modules/local/stats.nf'
include { MAKE_HEADERFILE       } from '../modules/local/header.nf'
include { MERGE_SORT ; REMOVE_DUPLICATES } from '../modules/local/fragments.nf'
include { MERGE_SORT_SAM        } from '../modules/local/bam.nf'

def buildFastqChannel() {
    if (!params.input) {
        exit(1, "No input samplesheet provided. Please specify with --input")
    }

    def readSuffix = "${params.readstr1 ?: ''}${params.ext ?: ''}"

    Channel.fromPath(params.input)
        .ifEmpty { exit(1, "Input samplesheet not found: ${params.input}") }
        .splitCsv(header: true)
        .map { row ->
            def firstEntry = row.values().find { it != null && it.toString().trim() }
            if (!firstEntry) {
                return null
            }
            if (firstEntry.toString().trim().startsWith('#')) {
                return null
            }
            return row
        }
        .filter { it != null }
        .ifEmpty { exit(1, "Samplesheet ${params.input} does not contain any usable records") }
        .map { row ->
            def sample = (row.sample ?: row.Sample ?: row.sample_id ?: row.Sample_ID)?.toString()?.trim()
            def name = (row.name ?: row.Name ?: row.library ?: row.Library)?.toString()?.trim()
            def fastq1 = (row.fastq_1 ?: row.fastq1 ?: row.read1 ?: row.r1)?.toString()?.trim()
            def fastq2 = (row.fastq_2 ?: row.fastq2 ?: row.read2 ?: row.r2)?.toString()?.trim()

            if (!sample) {
                exit(1, "Samplesheet is missing a 'sample' column entry")
            }
            if (!fastq1 || !fastq2) {
                exit(1, "Samplesheet entry for sample '${sample}' is missing FASTQ paths")
            }

            def fq1Path = file(fastq1)
            def fq2Path = file(fastq2)

            if (!fq1Path.exists()) {
                exit(1, "FASTQ file not found: ${fastq1}")
            }
            if (!fq2Path.exists()) {
                exit(1, "FASTQ file not found: ${fastq2}")
            }

            if (!name) {
                def fqName = fq1Path.getFileName().toString()
                if (readSuffix && fqName.endsWith(readSuffix)) {
                    name = fqName[0..<(fqName.length() - readSuffix.length())]
                }
                else {
                    name = fqName.replaceAll(/(_R1|_1)([^_\.]*)?$/, '')
                }
            }

            tuple(sample, name, fq1Path, fq2Path)
        }
        .share()
}

def validateParameters() {

    if (!params.readstr1) {
        params.readstr1 = '_R1'
    }

    if (!params.readstr2) {
        params.readstr2 = '_R2'
    }

    if (!params.reference) {
        exit(1, "Parameter --reference is required")
    }

    def referencePath = file(params.reference)
    if (!referencePath.exists()) {
        exit(1, "Reference FASTA not found: ${params.reference}")
    }

    params.reference = referencePath

    if (!params.genome_id) {
        exit(1, "Parameter --genome_id is required")
    }

    if (params.save_merged_nodups_bam == null) {
        params.save_merged_nodups_bam = true
    }
    else if (params.save_merged_nodups_bam instanceof Boolean) {
    }
    else if (params.save_merged_nodups_bam instanceof Number) {
        params.save_merged_nodups_bam = params.save_merged_nodups_bam as Integer != 0
    }
    else {
        def value = params.save_merged_nodups_bam.toString().toLowerCase()
        if (value in ['true', '1', 'yes', 'y']) {
            params.save_merged_nodups_bam = true
        }
        else if (value in ['false', '0', 'no', 'n']) {
            params.save_merged_nodups_bam = false
        }
        else {
            exit(1, "Parameter --save_merged_nodups_bam must be true or false")
        }
    }

    params.site = params.site.toString().toLowerCase()

    if (params.site == 'hindiii') {
        params.ligation = "AAGCTAGCTT"
    }
    else if (params.site == "dpnii") {
        params.ligation = "GATCGATC"
    }
    else if (params.site == 'mboi') {
        params.ligation = "GATCGATC"
    }
    else if (params.site == 'ncoi') {
        params.ligation = "CCATGCATGG"
    }
    else if (params.site == "arima") {
        params.ligation = "'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'"
    }
    else if (params.site == "none") {
        params.ligation = "XXXX"
    }
    else {
        exit(1, "Parameter --site must be one of: hindiii, mboi, dpnii, ncoi, arima, none")
    }

    if (!params.site_file) {
        exit(1, "Parameter --site_file is required")
    }
    else {
        if (!file(params.site_file).exists()) {
            exit(1, "Restriction site file not found: ${params.site_file}")
        }
    }

    if (!params.resolutions) {
        params.resolutions = '2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100'
    }

    if (!params.nofrag) {
        params.nofrag = true
    }
    else {
        def value = params.nofrag.toString().toLowerCase()
        if (value in ['true', '1', 'yes', 'y']) {
            params.nofrag = true
        }
        else if (value in ['false', '0', 'no', 'n']) {
            params.nofrag = false
        }
        else {
            exit(1, "Parameter --nofrag must be true or false")
        }
    }

    if (!params.justexact) {
        params.justexact = false
    }
    else {
        def value = params.justexact.toString().toLowerCase()
        if (value in ['true', '1', 'yes', 'y']) {
            params.justexact = true
        }
        else if (value in ['false', '0', 'no', 'n']) {
            params.justexact = false
        }
        else {
            exit(1, "Parameter --justexact must be true or false")
        }
    }

    if (!params.ext) {
        params.ext = '.fastq.gz'
    }
}

workflow NFCORE_JUICER {
    main:
    validateParameters()
    fastq_pairs = buildFastqChannel()

    out_header = MAKE_HEADERFILE(
        fastq_pairs.map { sample, _name, _read1, _read2 -> sample }.distinct()
    )

    frag_results = process_fragments(fastq_pairs)

    chimeric_reads = frag_results.chimeric_output
    sorted_fragments = frag_results.sorted_fragments

    sort_files_by_sample = sorted_fragments.groupTuple(by: 0)

    merged_sort = MERGE_SORT(sort_files_by_sample)

    merged_nodups = REMOVE_DUPLICATES(merged_sort)
    nodups = merged_nodups.map { sample, merged_nodups_txt, _dups_txt, _opt_dups_txt ->
        tuple(sample, merged_nodups_txt)
    }

    norm_sam_by_sample = chimeric_reads
        .map { sample, _name, _norm_txt, _abnorm_sam, _unmapped_sam, norm_sam, _norm_res_txt ->
            tuple(sample, norm_sam)
        }
        .groupTuple(by: 0)

    if (params.save_merged_nodups_bam) {
        merged_sam = MERGE_SORT_SAM(norm_sam_by_sample)
        def merged_nodups_by_sample = nodups.join(merged_sam)
        SAM_TO_BAM(REMOVE_DUPLICATES_SAM(merged_nodups_by_sample))
    }

    chimeric_by_sample = chimeric_reads
        .map { sample, _name, _norm_txt, abnorm_sam, unmapped_sam, _norm_sam, norm_res_txt ->
            tuple(sample, norm_res_txt, abnorm_sam, unmapped_sam)
        }
        .groupTuple(by: 0)

    stats_input = out_header
        .join(chimeric_by_sample)
        .join(nodups)

    stats_output = STATS(stats_input)

    hic_input = stats_output
        .map { sample, inter, inter_30, inter_hists, _collisions, _abnorm_sam, _unmapped_sam ->
            tuple(sample, inter, inter_30, inter_hists)
        }
        .join(nodups)

    hic_out_ch = GEN_HIC_FILES(hic_input)

    inter_30_hic = hic_out_ch.map { sample, _inter_hic, inter_30_hic_file, _inter_30_hists ->
        tuple(sample, inter_30_hic_file)
    }

    postprocessing(inter_30_hic)

    emit:
    stats = stats_output
    hic   = hic_out_ch
}
