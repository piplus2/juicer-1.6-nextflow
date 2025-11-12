# nf-core/juicer

An nf-core style port of the Juicer 1.6 Hi-C processing workflow. The pipeline orchestrates the complete Juicer analysis using Nextflow DSL2 modules and supports both CPU and GPU enabled stages such as read alignment, fragment processing, `.hic` generation, Arrowhead and HiCCUPS loop calling, and motif discovery.

## Key features

- **nf-core layout** – Standardised project structure with `main.nf`, `workflow/`, and modular DSL2 components for easier maintenance and customisation.
- **Samplesheet driven inputs** – Supply paired-end FASTQ files through a CSV samplesheet for reproducible and batch-friendly submissions.
- **Configurable profiles** – `standard` (local) and `hpc` (PBS Pro) profiles with consistent process labels and containers, plus an extendable `test` profile stub.
- **GPU-aware post-processing** – HiCCUPS, Arrowhead, and motif discovery steps honour GPU resource labels while falling back gracefully if CUDA is unavailable.

## Quick start

1. Install [Nextflow](https://www.nextflow.io/) (>= 23.10.0) and ensure Java 11+ is available.
2. Clone the repository and change into the project directory:

   ```bash
   git clone https://github.com/nf-core/juicer.git
   cd juicer-1.6-nextflow
   ```

3. Prepare a CSV samplesheet with the following columns:

   ```text
   sample,name,fastq_1,fastq_2
   sampleA,sampleA,/data/fastq/sampleA_R1.fastq.gz,/data/fastq/sampleA_R2.fastq.gz
   sampleB,sampleB,/data/fastq/sampleB_R1.fastq.gz,/data/fastq/sampleB_R2.fastq.gz
   ```

4. Launch the pipeline:

   ```bash
   nextflow run . \
     --input samplesheet.csv \
     --reference /path/to/genome.fa \
     --genome_id mm10 \
     --outdir results \
     --site arima \
     --motif_dir /path/to/motif_directory \
     -profile standard
   ```

   Use `-profile hpc` to activate the bundled PBS Pro configuration or extend the configuration files in `conf/` for your environment.

## Parameters and configuration

All default parameters are defined in `nextflow.config`. Important options include:

- `--input` – CSV samplesheet describing paired FASTQ files.
- `--reference` – Path (local or remote) to the reference genome FASTA. The file is staged by Nextflow before execution.
- `--genome_id` – Genome identifier passed to Juicer tools.
- `--site` / `--site_file` – Restriction enzyme label or explicit restriction site file. The pipeline resolves known motifs automatically when possible.
- `--motif_dir` – Directory containing motif files for loop annotation (optional).
- `--outdir` – Destination for all published results (default: `./results`).
- `--threads` – Default CPU allocation for processes carrying the `hpc` label.
- `--java_mem` – Maximum Java heap size supplied to Juicer tools.
- `--align_cpus` / `--align_memory_gb` / `--align_time` – Fine-tune resources for `BWA_ALIGN` (defaults: `min(--threads,16)` CPUs, `64 GB`, `48h`).
- `--default_cpus` / `--default_memory_gb` / `--default_time` – Baseline resources for non-GPU processes without explicit overrides (defaults: `4`, `16 GB`, `48h`).
- `--default_gpu_cpus` / `--default_gpu_memory_gb` / `--default_gpu_time` / `--gpu_cluster_options` – Baseline resources for GPU-labelled stages (defaults: `8`, `64 GB`, `48h`, and PBS `-l select=1:ncpus=24:ngpus=1`).
- Per-process overrides follow the pattern `--<process>_cpus`, `--<process>_memory_gb`, and `--<process>_time` (e.g. `--stats_memory_gb 48` or `--sam_to_bam_cpus 12`). See `conf/process_resources.config` for the full list of keys.
- `--save_merged_nodups_bam` – Toggle creation of the deduplicated `merged_nodups.bam` file (default: enabled).

Profiles are defined in `nextflow.config`:

- `standard` – Local execution with default resource requests.
- `hpc` – PBS Pro executor configuration mirroring the original cluster deployment.
- `test` – Stub profile pointing to `assets/test_samplesheet.csv`; replace the example paths with real small test data before running.

## Pipeline overview

The DSL2 workflow mirrors the original Juicer stages:

1. **FASTQ validation and ligation counting** – Ensures mate pairs exist, counts ligation motifs, and prepares fragment inputs.
2. **Alignment and fragment processing** – Aligns reads with `bwa-mem2`, generates fragment files, and sorts them by genome coordinates.
3. **Duplicate removal and BAM conversion** – Merges fragment outputs, removes duplicates, and converts filtered SAM files to BAM using `samtools`.
4. **Statistics and `.hic` production** – Generates standard Juicer statistics (`inter.txt`, `inter_30.txt`, etc.) and produces both standard and 30 kb `.hic` files.
5. **Post-processing (GPU-enabled)** – Runs Arrowhead, HiCCUPS, and optional motif discovery against provided motif libraries.

Outputs are organised by sample inside `--outdir`, mirroring the legacy Juicer folder structure (`aligned/`, `splits/`, and post-processing directories).

## Resource heuristics

`conf/process_resources.config` captures the per-process resource matrix used by every profile. CPUs and memory are derived from lightweight heuristics bounded by the user-facing parameters listed above. The defaults are summarised below:

| Process | Default CPUs | Default Memory | Notes |
|---------|--------------|----------------|-------|
| `COUNT_LIGATIONS` | 2 | 4 GB | Ligations counting through `paste`/`grep`. |
| `BWA_ALIGN` | `--align_cpus` (defaults to `min(--threads, 16)`) | 64 GB (override via `--align_memory_gb`) | Heavy bwa-mem2 alignment. |
| `FRAGMENT` | 2 | 8 GB | `fragment.pl` or AWK based fragment builder. |
| `SORT` | 4 | 32 GB | Large coordinate sort (`sort -m`). |
| `CHIMERIC` | 4 | 16 GB | Juicer `chimeric_blacklist.awk` filtering. |
| `MERGE_SORT` | 4 | 16 GB | Multi-way merge of fragment files. |
| `REMOVE_DUPLICATES` | 4 | 12 GB | `dups.awk` duplicate culling. |
| `REMOVE_DUPLICATES_SAM` | 4 | 12 GB | Read-name based SAM filtration. |
| `SAM_TO_BAM` | 8 | 24 GB | `samtools view` conversion (`-@ cpus`). |
| `MERGE_SORT_SAM` | 8 | 32 GB | `samtools merge|sort` of SAM chunks. |
| `MAKE_HEADERFILE` | 1 | 2 GB | Metadata banner writer. |
| `STATS` | 4 | 32 GB | Juicer statistics (`juicer_tools`) plus AWK summaries. |
| `GEN_HIC_FILES` | 6 | 48 GB | Two passes of `juicer_tools pre` and `statistics.pl`. |
| `JUICER_TOOLS_PRE_Q1/Q30` | 4 | 32 GB | Standalone `.hic` creation. |
| `HIC_STATS` | 4 | 24 GB | Generates `inter_30_hists.m`. |
| `ARROWHEAD` | `--default_gpu_cpus` (default 8) | `--default_gpu_memory_gb` (default 64 GB) | GPU-friendly domain calling. |
| `HICCUPS` | `max(8, --default_gpu_cpus)` | 96 GB (override via `--hiccups_memory_gb`) | GPU loop calling. |
| `MOTIF_FINDER` | `--default_gpu_cpus` (default 8) | 64 GB (override via `--motif_finder_memory_gb`) | `juicer_tools apa` + motif search. |

Any process not explicitly listed falls back to `--default_cpus` / `--default_memory_gb` / `--default_time`. Override individual entries directly from the CLI, for example:

```bash
nextflow run . \
  --sam_to_bam_cpus 12 \
  --sam_to_bam_memory_gb 48 \
  --hiccups_memory_gb 128
```

## Additional documentation

Extended usage notes and parameter descriptions are available in `docs/usage.md`. Update or extend the documentation as you adapt the pipeline to your environment.

## Acknowledgements

- Original Nextflow implementation by **Paolo Inglese** (Fondazione Istituto Italiano di Tecnologia).
- Ported to the nf-core style layout by the OpenAI Autopilot.
- Juicer tools are distributed by the Aiden Lab.

## License

Released under the MIT License. See `LICENSE` for full details.
