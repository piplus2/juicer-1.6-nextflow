# Usage guide

This document supplements the main `README.md` with additional tips for operating the nf-core flavoured Juicer pipeline.

## Samplesheet format

The workflow expects a comma-separated file with at least the columns below. Column names are case-insensitive and additional columns are ignored.

| Column    | Description                                                                 |
|-----------|-----------------------------------------------------------------------------|
| `sample`  | Biological sample identifier.                                               |
| `name`    | Library or lane identifier used to tag intermediate files.                  |
| `fastq_1` | Absolute or relative path to the R1 FASTQ file.                              |
| `fastq_2` | Absolute or relative path to the R2 FASTQ file.                              |

If the `name` column is omitted the pipeline derives it from the R1 file name by removing `--readstr1` and `--ext` (defaults `_R1_001` and `.fastq.gz`).

Example:

```text
sample,name,fastq_1,fastq_2
patient1,patient1_lane1,data/patient1_L001_R1.fastq.gz,data/patient1_L001_R2.fastq.gz
patient1,patient1_lane2,data/patient1_L002_R1.fastq.gz,data/patient1_L002_R2.fastq.gz
```

## Common parameters

| Parameter | Purpose |
|-----------|---------|
| `--reference` | Reference genome FASTA. The path is validated before any process starts. |
| `--genome_id` | Genome label passed to Juicer tools (e.g. `hg38`, `mm10`). |
| `--site` | Restriction enzyme name (e.g. `arima`, `hindiii`). Automatically sets ligation motifs. |
| `--site_file` | Override restriction site file if the automatic lookup does not match your genome. |
| `--motif_dir` | Directory with motif definitions for the optional motif finding step. |
| `--nofrag` | Skip fragment delimitation (`1` by default). Set to `0` to enable fragment-based maps. |
| `--resolutions` | Comma-separated list of map resolutions passed to HiCCUPS and `.hic` generation. |
| `--save_merged_nodups_bam` | Set to 'false' to skip generating the deduplicated `merged_nodups.bam` output (default: 'true'). |

## Profiles

- `standard`: Local execution suitable for development and small runs.
- `hpc`: PBS Pro profile that mirrors the IIT deployment. Adjust queue names, scratch behaviour, and `clusterOptions` as needed.
- `test`: Convenience profile pointing to `assets/test_samplesheet.csv`. Replace placeholder paths before running.

## Resource tuning

Resource requests are derived from the matrix in `conf/process_resources.config`. Customize them without editing the file by passing CLI flags:

- `--default_cpus`, `--default_memory_gb`, `--default_time`: Baseline for all CPU-only processes.
- `--default_gpu_cpus`, `--default_gpu_memory_gb`, `--default_gpu_time`, `--gpu_cluster_options`: Baseline for GPU-labelled stages.
- Per-process overrides follow the convention `--<process>_cpus`, `--<process>_memory_gb`, `--<process>_time` (e.g. `--merge_sort_sam_cpus 12`, `--gen_hic_files_memory_gb 56`).
- Alignment-specific helpers: `--align_cpus`, `--align_memory_gb`, `--align_time`.

Example CLI overrides:

```bash
nextflow run . \
  --align_cpus 24 \
  --align_memory_gb 96 \
  --stats_memory_gb 48 \
  --hiccups_memory_gb 128
```

## Containers and dependencies

The pipeline assumes Juicer helper scripts and binaries are available in `bin/` (added to the `PATH` via `conf/base.config`). Containers are only configured for `samtools`; feel free to extend the configuration with Singularity/Docker images for the remaining tools.

## GPU notes

HiCCUPS and motif discovery require GPUs. The `gpu` process label requests a GPU resource on the HPC profile. When CUDA is unavailable the pipeline exits with a clear error message.

## Extending the pipeline

- Add new modules under `modules/local/` or import nf-core/community modules using the standard layout.
- Additional subworkflows can live under `subworkflows/local/` and be included from `workflow/main.nf`.
- Custom profiles should be added to the `profiles` block in `nextflow.config` and optionally to separate files in `conf/`.

Happy processing!
