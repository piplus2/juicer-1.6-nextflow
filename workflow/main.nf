nextflow.enable.dsl = 2

include { process_fragments     } from '../subworkflows/local/align_reads'
include { postprocessing        } from '../subworkflows/local/postproc'
include { REMOVE_DUPLICATES_SAM } from '../modules/local/bam/remove_dups_sam'
include { SAM_TO_BAM            } from '../modules/local/bam/convert_to_bam'
include { hic                   } from '../subworkflows/local/gen_hic'
include { STATS                 } from '../modules/local/stats'
include { MAKE_HEADERFILE       } from '../modules/local/header'
include { MERGE_SORT            } from '../modules/local/merge_sort'
include { REMOVE_DUPLICATES     } from '../modules/local/remove_dups'
include { MERGE_SORT_SAM        } from '../modules/local/bam/merge_and_sort'

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
}

def validateParameters() {
    if (!params.outdir) {
        params.outdir = "${launchDir}/results"
    }

    if (!params.output_dir) {
        params.output_dir = params.outdir
    }

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
    else {
        params.genomeID = params.genome_id
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

    if (params.use_gpu == null) {
        params.use_gpu = true
    }
    else if (params.use_gpu instanceof Boolean) {
    }
    else if (params.use_gpu instanceof Number) {
        params.use_gpu = params.use_gpu as Integer != 0
    }
    else {
        def value = params.use_gpu.toString().toLowerCase()
        if (value in ['true', '1', 'yes', 'y']) {
            params.use_gpu = true
        }
        else if (value in ['false', '0', 'no', 'n']) {
            params.use_gpu = false
        }
        else {
            exit(1, "Parameter --use_gpu must be true or false")
        }
    }

    if (!params.site) {
        exit(1, "Parameter --site must be specified")
    }

    def validSites = ['hindiii', 'mboi', 'dpnii', 'ncoi', 'arima', 'none']
    if (!(params.site in validSites)) {
        exit(1, "Parameter --site must be one of: ${validSites.join(', ')}")
    }

    if (!params.ligation) {
        exit(1, "Parameter --ligation is required (normally derived automatically from --site)")
    }

    if (!params.site_file) {
        exit(1, "Parameter --site_file is required")
    }
    else {
        if (!file(params.site_file).exists()) {
            exit(1, "Restriction site file not found: ${params.site_file}")
        }
    }

    if (!params.resolutions.matches(/^(\d+,)*\d+$/)) {
        exit(1, "Parameter --resolutions must be a comma-separated list of integers")
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

    if (params.skip_motif_finder == null) {
        params.skip_motif_finder = false
    }
    else if (params.skip_motif_finder instanceof Boolean) {
    }
    else if (params.skip_motif_finder instanceof Number) {
        params.skip_motif_finder = params.skip_motif_finder as Integer != 0
    }
    else {
        def value = params.skip_motif_finder.toString().toLowerCase()
        if (value in ['true', '1', 'yes', 'y']) {
            params.skip_motif_finder = true
        }
        else if (value in ['false', '0', 'no', 'n']) {
            params.skip_motif_finder = false
        }
        else {
            exit(1, "Parameter --skip_motif_finder must be true or false")
        }
    }
}


def displayInfo() {
    println("================ Juicer Pipeline ================")
    println("Reference genome:       ${params.reference}")
    println("Genome ID:              ${params.genome_id}")
    println("Restriction site:       ${params.site} (${params.ligation})")
    println("Restriction site file:  ${params.site_file}")
    println("Output directory:       ${params.outdir}")
    println("Save merged no-dup BAM: ${params.save_merged_nodups_bam}")
    println("Read 1 suffix:         '${params.readstr1}'")
    println("Read 2 suffix:         '${params.readstr2}'")
    println("FASTQ extension:        '${params.ext}'")
    println("No fragment filtering:  ${params.nofrag}")
    println("Just exact matches:     ${params.justexact}")
    println("Resolutions:            ${params.resolutions}")
    println("================================================")
}


workflow NFCORE_JUICER {
    main:
    validateParameters()
    displayInfo()
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

    if (params.save_merged_nodups_bam) {
        norm_sam_by_sample = chimeric_reads
            .map { sample, _name, _norm_txt, _abnorm_sam, _unmapped_sam, norm_sam, _norm_res_txt ->
                tuple(sample, norm_sam)
            }
            .groupTuple(by: 0)

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
        .join(merged_nodups)

    site_file = file(params.site_file)
    stats_output = STATS(stats_input, site_file)

    hic_input = stats_output
        .map { sample, inter_txt, inter_30_txt, inter_hists_m, _collisions, _abnorm_sam, _unmapped_sam ->
            tuple(sample, inter_txt, inter_30_txt, inter_hists_m)
        }
        .join(nodups)

    hic_out_ch = hic(hic_input)

    postprocessing(hic_out_ch)

    emit:
    stats = stats_output
    hic   = hic_out_ch
}
