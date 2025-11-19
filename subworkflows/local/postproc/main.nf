include { ARROWHEAD    } from '../../../modules/local/juicer_tools/arrowhead'
include { HICCUPS      } from '../../../modules/local/juicer_tools/hiccups'
include { MOTIF_FINDER } from '../../../modules/local/juicer_tools/motif_finder'


workflow postprocessing {
    take:
    inter_30_hic

    main:

    ARROWHEAD(inter_30_hic)
    loops_dir = HICCUPS(inter_30_hic)
    // Keep the loops dirs that contain the file merged_loops.bedpe
    good_loops_dir = loops_dir.filter { _sample, _inter_30, _loops_dir ->
        file("${_loops_dir}/merged_loops.bedpe").exists()
    }

    if (params.skip_motif_finder == false) {
        motif_dir = Channel.fromPath(params.motif_dir, checkIfExists: true).first()
        MOTIF_FINDER(good_loops_dir, motif_dir)
    }

    emit:
    null
}
