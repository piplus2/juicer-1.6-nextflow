# nf-core/juicer

An nf-core inspired Nextflow DSL2 rewrite of the Juicer 1.6 Hi-C processing system. The project wraps the complete Juicer workflow—from FASTQ validation to Arrowhead, HiCCUPS, and motif discovery—inside a modular Nextflow pipeline so that large Hi-C projects can be reproduced on laptops, workstations, or schedulers with consistent configuration.

## What this repository contains

- `main.nf` and `workflow/` – Nextflow DSL2 entry point and workflow graph.
- `modules/` and `subworkflows/` – Reusable tasks mirroring every Juicer stage (alignment, fragment processing, `.hic` generation, post-processing, etc.).
- `conf/` – Base resources plus ready-to-use profiles for local (`standard`), PBS Pro HPC (`hpc`), tests, and optional `conda` / `singularity` runtimes.
- `assets/` – Example data such as `test_samplesheet.csv`.
- `bin/` – Juicer helper scripts (AWK, Perl, bash) referenced by the modules.
- `docs/usage.md` – Extended usage tips and parameter explanations.

## Requirements

- [Nextflow](https://www.nextflow.io/) `>= 23.10.0`
- Java 11+
- Juicer helper scripts (already vendored in `bin/`)
- `bwa-mem2`, `samtools`, `juicer_tools`, GNU coreutils, Perl, and AWK
- Optional CUDA-capable GPUs for HiCCUPS and motif discovery stages

Install tools through modules, Conda, Singularity, system packages, or your cluster environment. The bundled profiles assume binaries are on `PATH`; adapt `conf/base.config` if you rely on containers or modules.

## Inputs

The pipeline accepts a comma-separated samplesheet (case-insensitive header names):

| Column    | Description                                                                 |
| --------- | --------------------------------------------------------------------------- |
| `sample`  | Biological sample identifier; controls the final directory structure.       |
| `name`    | Technical replicate / lane identifier used when tagging intermediate files. |
| `fastq_1` | Absolute path to the R1 FASTQ file.                                         |
| `fastq_2` | Absolute path to the R2 FASTQ file.                                         |

Example (`assets/test_samplesheet.csv` mirrors this layout):

```text
sample,name,fastq_1,fastq_2
sampleA,sampleA_lane1,/data/hi-c/sampleA_L001_R1.fastq.gz,/data/hi-c/sampleA_L001_R2.fastq.gz
sampleA,sampleA_lane2,/data/hi-c/sampleA_L002_R1.fastq.gz,/data/hi-c/sampleA_L002_R2.fastq.gz
```

## Quick start

1. Install Nextflow and Java, then clone the repository:

   ```bash
   git clone git@github.com:piplus2/juicer-1.6-nextflow.git
   cd juicer-1.6-nextflow
   ```

2. Edit a samplesheet as described above.
3. Launch the workflow:

   ```bash
   nextflow run . \
     --input samplesheet.csv \
     --reference /path/to/genome.fa \
     --genome_id hg38 \
     --site arima \
     --motif_dir /path/to/motifs \
     --outdir results \
     -profile standard
   ```

Add `-profile hpc` for the PBS Pro executor, combine `conda` or `singularity` with other profiles (e.g. `-profile standard,singularity`), or create custom configs in `conf/`.

## Pipeline overview

1. **FASTQ validation & ligation counting** – Ensures paired files exist, counts ligation motifs, and prepares fragment metadata.
2. **Alignment & fragment building** – Runs `bwa-mem2` against the supplied reference, generates sorted fragment files, and handles split/ambiguous reads.
3. **Merging, duplicate removal & BAM creation** – Performs juicer-style merges, deduplication, and exports filtered BAMs (`merged_nodups.bam` optional).
4. **Statistics & `.hic` generation** – Produces Juicer statistics (`inter.txt`, `inter_30.txt`, histograms) and `.hic` files at requested resolutions.
5. **Post-processing & motif discovery** – Executes Arrowhead, HiCCUPS, and optional motif scanning using GPU resources when available.

Results mimic the classic Juicer layout under `--outdir` (e.g. `aligned/`, `splits/`, `juicer_tools/arrowhead`, `juicer_tools/hiccups`, `juicer_tools/motifs`).

## Configuration essentials

All default parameters live in `nextflow.config`. Frequently tuned options:

| Parameter                  | Purpose                                                             |
| -------------------------- | ------------------------------------------------------------------- |
| `--input`                  | Samplesheet described above.                                        |
| `--reference`              | Reference genome FASTA (local path or remote URL).                  |
| `--genome_id`              | Genome identifier passed to Juicer (hg38, mm10, etc.).              |
| `--site` / `--site_file`   | Restriction enzyme shorthand or explicit restriction site BED.      |
| `--motif_dir`              | Repository of motif definitions for the optional motif finder.      |
| `--nofrag`                 | Set to `0` to enable fragment-based maps (default `1`).             |
| `--resolutions`            | Comma-separated list of map resolutions used by `.hic` and HiCCUPS. |
| `--threads`                | Baseline CPU count that informs resource heuristics (default `16`). |
| `--java_mem`               | Max heap supplied to `juicer_tools`.                                |
| `--save_merged_nodups_bam` | Toggle writing `merged_nodups.bam` (default `true`).                |

### Profiles

- `standard` – Local execution with conservative CPU/memory defaults.
- `hpc` – PBS Pro template (swap queues and directives for SLURM or other schedulers).
- `conda` / `singularity` – Enable each runtime; combine with other profiles as needed.

### Resource overrides

Resource classes (`smallcpu`, `mediumcpu`, `highcpu`, `gpu`) are defined in `conf/process_resources.config`. You can override them without editing config files:

```bash
nextflow run . \
  --highcpu_cpus 32 \
  --highcpu_memory_gb 128 \
  --align_cpus 24 \
  --hiccups_memory_gb 128
```

Per-process overrides follow the `--<process>_cpus|memory_gb|time` convention (for example `--sam_to_bam_cpus 12`). See `docs/usage.md` for an expanded matrix and helper flags such as `--align_memory_gb`.

## Containers and dependencies

The repository includes required Juicer helper scripts in the `bin/` directory. By default, the pipeline uses containers for `bwa-mem2`, samtools`and`juicer_tools`. To run additional tools inside Singularity or Docker, extend the definitions in `conf/containers.config`. For cluster environments that use environment modules, set the appropriate `beforeScript`or`module` directives in your profile.

## Testing and development

- Run smoke tests with `nextflow run . -profile test`.
- Import additional DSL2 modules under `modules/local/` or plug in nf-core/community modules for new features.
- Add bespoke profiles inside `nextflow.config` or to dedicated files under `conf/`.

## Documentation

`docs/usage.md` contains supplementary guidance covering samplesheet nuances, CLI parameter tips, and GPU resource notes. Update both this README and the usage guide when you add new features so users can see what has changed at a glance.

## Acknowledgements

- Original Juicer code and `juicer_tools` were created by Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, Erez Lieberman Aiden, and colleagues at the Aiden Lab (see **Durand et al., Cell Systems 2016** for the canonical citation). Please cite their paper when publishing results produced with this workflow.
- This Nextflow port was developed by **Paolo Inglese**.
- Portions of the implementation borrow directly from the original Juicer helper scripts distributed under the Juicer license; those files retain their original headers.

## License

Released under the MIT License. See `LICENSE` for the full text.
