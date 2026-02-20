# Nextflow Tutorial

A hands-on tutorial for learning Nextflow, a workflow management system for bioinformatics pipelines. This tutorial progresses from basic "Hello World" examples to a complete parallel BLAST workflow.

## Table of Contents
- [Setup](#setup)
- [Tutorial Overview](#tutorial-overview)
- [Tutorial Examples](#tutorial-examples)
- [Running the Examples](#running-the-examples)
- [Project Structure](#project-structure)
- [Resources](#resources)

## Setup

### Prerequisites
- Access to an HPC environment with Nextflow and BLAST+ installed
- Basic knowledge of command-line operations
- Understanding of FASTA format and BLAST basics (for later examples)

### Initial Configuration

1. **Load Nextflow module:**
   ```bash
   module load nextflow
   ```

2. **Set Nextflow home directory:**
   ```bash
   export NXF_HOME='/work/idoerg/<netid>/.nextflow'
   ```
   Replace `<netid>` with your actual network ID.

3. **Verify installation:**
   ```bash
   nextflow -h
   ```

## Tutorial Overview

This tutorial demonstrates key Nextflow concepts through seven progressively complex examples:

1. **Basic output** - Writing to stdout
2. **File output** - Creating output files
3. **Publishing results** - Using `publishDir` directive
4. **Multi-language support** - Running Python scripts
5. **Multiple processes** - Chaining multiple processes
6. **Real-world pipeline** - Complete BLAST workflow
7. **Parallel processing** - Splitting work across chunks

## Tutorial Examples

### 01_hello.nf - Basic Hello World
**Concepts:** Process definition, stdout, basic workflow

A simple introduction to Nextflow that prints a message to standard output.

```bash
nextflow run code/01_hello.nf
```

**Key concepts:**
- `nextflow.enable.dsl=2` - Enables DSL2 syntax
- `process` blocks define computational tasks
- `output` directive specifies what the process produces
- `workflow` block orchestrates process execution
- `.view()` displays output to console

---

### 02_hellowrite2file.nf - Writing to Files
**Concepts:** File outputs, output paths

Demonstrates how to write process output to a file instead of stdout.

```bash
nextflow run code/02_hellowrite2file.nf
```

**Key concepts:**
- `path 'result.txt'` specifies file output
- Output files are stored in Nextflow's `work/` directory by default

---

### 03_hellofile.nf - Publishing Results
**Concepts:** publishDir, output management

Shows how to copy results to a designated output directory.

```bash
nextflow run code/03_hellofile.nf
```

**Key concepts:**
- `publishDir` directive copies results to specified location
- `mode: 'copy'` creates a copy (alternatives: 'symlink', 'move')
- Results appear in the `output/` directory

---

### 04_hellopython.nf - Multi-Language Scripts
**Concepts:** Shebang, Python in Nextflow

Demonstrates running Python code within a Nextflow process.

```bash
nextflow run code/04_hellopython.nf
```

**Key concepts:**
- `#!/usr/bin/python` shebang specifies interpreter
- Any scripting language can be used (Python, R, Perl, etc.)
- Script content goes in the triple-quoted string

---

### 05_hellomultiprocess.nf - Multiple Processes
**Concepts:** Process composition, parallel execution

Shows how to define and run multiple independent processes.

```bash
nextflow run code/05_hellomultiprocess.nf
```

**Key concepts:**
- Multiple processes can be defined in one workflow
- Processes without dependencies run in parallel
- Each process has its own execution environment

---

### 06_blast.nf - Complete BLAST Pipeline
**Concepts:** Input parameters, process chaining, data dependencies

A real bioinformatics pipeline that downloads a database, creates BLAST indices, and performs sequence alignment.

```bash
nextflow run code/06_blast.nf -profile hpc_modules
```

**Key concepts:**
- `params` define configurable workflow parameters
- `input` and `output` directives chain processes together
- Processes execute only when their inputs are ready
- Profile configuration (in `nextflow.config`) manages environment

**Pipeline steps:**
1. `DOWNLOAD_UNIPROT_FASTA` - Downloads yeast proteins from UniProt
2. `MAKE_BLAST_DB` - Creates BLAST database from downloaded sequences
3. `BLASTP` - Runs protein BLAST search with query sequences

---

### 07_blastparallel.nf - Parallel BLAST
**Concepts:** Channel operations, data parallelization, result aggregation

Extends the BLAST pipeline with parallel processing by splitting the query file into chunks.

```bash
nextflow run code/07_blastparallel.nf -profile hpc_modules
```

**Key concepts:**
- `Channel.fromPath()` creates channels from files
- `.splitFasta()` divides FASTA files into chunks
- `.map()` transforms channel data
- `.collect()` aggregates results from parallel processes
- `tuple` passes multiple values together

**Pipeline steps:**
1. `DOWNLOAD_UNIPROT_FASTA` - Downloads database
2. `MAKE_BLAST_DB` - Creates BLAST indices
3. `BLASTP` - Runs BLAST on each chunk in parallel
4. `MERGE_BLAST` - Combines all chunk results into single file

**Performance benefit:** Processing 5 sequences at once (configurable via `params.chunkSize`) allows parallel execution and faster completion.

## Running the Examples

### Basic Execution
```bash
# Run without profile (for simple examples)
nextflow run code/01_hello.nf

# Run with HPC modules profile (for BLAST examples)
nextflow run code/06_blast.nf -profile hpc_modules
```

### Common Options
```bash
# Resume from last successful step
nextflow run code/07_blastparallel.nf -resume -profile hpc_modules

# Override parameters
nextflow run code/07_blastparallel.nf --chunkSize 10 -profile hpc_modules

# View execution timeline
nextflow run code/07_blastparallel.nf -profile hpc_modules -with-timeline timeline.html
```

### Cleaning Up
```bash
# Remove work directory after successful run
rm -rf work/

# Clean up Nextflow cache
nextflow clean -f
```

### Performance Comparison: Nextflow vs Bash

To compare execution times between Nextflow and traditional bash scripts, you can use the `time` command:

#### Timing the Bash Script (06_blast.sh)
```bash

# Time the bash script execution
time ./code/06_blast.sh
```

#### Timing the Nextflow Pipeline (06_blast.nf)
```bash
# Time the Nextflow pipeline
time nextflow run code/06_blast.nf -profile hpc_modules

# For subsequent runs (using cache)
time nextflow run code/06_blast.nf -profile hpc_modules -resume
```

#### Timing Output Explained
The `time` command shows three values:
- **real** - Total wall clock time (actual elapsed time)
- **user** - CPU time spent in user mode
- **sys** - CPU time spent in kernel mode

#### Example Comparison
```bash
# Clean run (no cached results)
rm -rf work/ work_bash/ results/

# Run bash version
echo "=== Bash Script ==="
time ./code/06_blast.sh

# Run Nextflow version
echo "=== Nextflow Pipeline ==="
time nextflow run code/06_blast.nf -profile hpc_modules
```

**Key Observations:**
- **First run**: Nextflow has overhead for workflow management
- **Cached runs**: Nextflow's `-resume` skips completed steps (huge time saver!)
- **Parallel workflows** (07_blastparallel.nf): Nextflow shows real advantage with parallelization
- **Reproducibility**: Nextflow tracks all intermediate steps automatically

#### Timing Parallel Execution
```bash
# Compare sequential vs parallel BLAST
echo "=== Sequential (06_blast.nf) ==="
time nextflow run code/06_blast.nf -profile hpc_modules

echo "=== Parallel (07_blastparallel.nf) ==="
time nextflow run code/07_blastparallel.nf -profile hpc_modules --chunkSize 5
```

## Project Structure

```
nextflow_tutorial/
├── README.md                    # This file
├── nextflow.config              # Configuration profiles
├── code/                        # Tutorial scripts
│   ├── 01_hello.nf             # Basic output
│   ├── 02_hellowrite2file.nf   # File output
│   ├── 03_hellofile.nf         # Publishing results
│   ├── 04_hellopython.nf       # Python integration
│   ├── 05_hellomultiprocess.nf # Multiple processes
│   ├── 06_blast.nf             # BLAST pipeline (Nextflow)
│   ├── 06_blast.sh             # BLAST pipeline (bash - for comparison)
│   └── 07_blastparallel.nf     # Parallel BLAST
├── data/                        # Input data
│   └── query.fasta             # Query sequences
├── output/                      # Simple example outputs
├── results/                     # BLAST results
├── work/                        # Nextflow working directory
└── work_bash/                   # Bash script working directory
```

## Configuration Profiles

The `nextflow.config` file defines two profiles:

### hpc_modules
For HPC environments using environment modules:
```bash
nextflow run code/06_blast.nf -profile hpc_modules
```
- Loads BLAST+ module before each task
- Uses local executor

### local
For local execution without modules:
```bash
nextflow run code/01_hello.nf -profile local
```
- Basic local execution
- Assumes tools are in PATH

## Key Nextflow Concepts

### Processes
- Self-contained computational tasks
- Define inputs, outputs, and script to execute
- Can use any scripting language

### Channels
- Asynchronous queues that connect processes
- Pass data between processes
- Support operations like map, filter, split, collect

### Workflows
- Define the execution logic
- Connect processes through channels
- Can be nested and modularized

### Directives
- `publishDir` - Where to save results
- `input` - What data the process receives
- `output` - What data the process produces
- `script` - The command(s) to execute

## Resources

### Official Documentation
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
- [nf-core](https://nf-co.re/) - Community curated pipelines

### Learning Resources
- [Nextflow Training](https://training.nextflow.io/)
- [Seqera Labs YouTube](https://www.youtube.com/@seqeralabs)

### Community
- [Nextflow Slack](https://www.nextflow.io/slack-invite.html)
- [nf-core Slack](https://nf-co.re/join/slack)

## Troubleshooting

### Common Issues

**Issue:** "command not found" errors
- **Solution:** Ensure modules are loaded or tools are in PATH. Use `-profile hpc_modules` for BLAST examples.

**Issue:** "work directory too large"
- **Solution:** Run `nextflow clean -f` to remove cached work files.

**Issue:** Pipeline fails partway through
- **Solution:** Use `-resume` flag to continue from last successful step.

**Issue:** Permission denied errors
- **Solution:** Check file permissions and ensure NXF_HOME is writable.

## Tips for Success

1. **Always use `-resume`** - Saves time by skipping completed steps
2. **Start simple** - Master basic examples before complex pipelines
3. **Check work directory** - Inspect intermediate files in `work/` for debugging
4. **Use `.view()`** - Add to channels to see what data is flowing through
5. **Read error messages** - Nextflow provides detailed logs and error traces

## Next Steps

After completing this tutorial:
1. Explore [nf-core pipelines](https://nf-co.re/pipelines) for production-ready workflows
2. Learn about [containers](https://www.nextflow.io/docs/latest/container.html) (Docker/Singularity) for reproducibility
3. Study [executor configuration](https://www.nextflow.io/docs/latest/executor.html) for HPC clusters
4. Build your own pipeline for your research needs

---
