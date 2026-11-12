# AbRaCD (AntiBody Random Coordinate Descent)

## Overview

AbRaCD (Antibody Random Coordinate Descent) is a sampling-based method for predicting antibody complementarity-determining region (CDR) loop structures. It generates and ranks multiple candidate conformations using antibody-specific structural constraints, then refines the selected models. The [AbRaCD web server](https://abracd.iqf.csic.es) provides access to the method.

This repository contains the scripts, data, and instructions needed to run AbRaCD locally and reproduce the results reported in the associated paper. It is intended to make the evaluation workflow straightforward to run and inspect.

> Martín-Hernández, I., Chacón, P., & López-Blanco, J. R. (manuscript in preparation). *AbRaCD: Fast sampling-based prediction of antibody CDR loops*. Bioinformatics Advances. [Publication details to be added when final.]

## Sections

- [Overview](#overview)
- [Repository directory guide](#repository-directory-guide)
- [Requirements](#requirements)
- [Installation](#installation)
  - [General Bash user without Conda](#general-bash-user-without-conda)
  - [RCD](#rcd)
- [Quick start](#quick-start)
- [Fast commands for reproducing benchmark paper results](#fast-commands-for-reproducing-benchmark-paper-results)
  - [AbRaCD Leakege Free Benchmark (ALFB)](#abracd-leakege-free-benchmark-alfb)
  - [RosettaBenchmark](#rosettabenchmark)
- [Repository contents](#repository-contents)
- [Troubleshooting](#troubleshooting)

---

- [Manual](AbRaCD_Manual.md)

## Repository directory guide

The main project directories are organized as follows. Paths are relative to the repository root.

```text
data/
├── Benchmarks/
│   ├── ALFB/
│   ├── RosettaBenchmark/
│   └── Therapeutics/
└── Maps/
paper_results/
scripts/
├── perl/
└── python/
src/
├── perl/
└── python/
```

- `data/`
  - `data/Benchmarks/` groups the available evaluation datasets (AFLB and Rosetta). Each dataset may include source Chothia structures, processed/renumbered structures with the RCD input lists, FASTA sequences, AlphaFold 3 model structures, and hybrid AF3–AbRaCD structures.
  - `data/Maps/` contains the binary maps required by the workflow. CDR is included and KORP, and Dunbrack should be copied after RCD instalation.
- `paper_results/` paper-related results.
- `scripts/` contains the command-line workflow scripts.
- `src/` contains supporting utilities used by the workflow.

## Requirements

- **Python 3.10 or later**, with `mpi4py`
- **Perl 5.32 or later**
- **MPI**, recommended
- **PyRosetta**, required for structure cleaning and loop refinement. [Download PyRosetta](https://www.pyrosetta.org/downloads)

Tested PyRosetta releases:

- `PyRosetta4.Release.python310.ubuntu.release-347`
- `PyRosetta4.Release.python310.ubuntu 2025.50 (weekly)`
- `PyRosetta4.Debug.python312.linux 2024.39 (weekly)`
- `PyRosetta4.Ubuntu.python314.Release 2026.39 (weekly)`

Other PyRosetta versions may work, but have not been tested. Compatibility is more likely with commonly used versions; very old versions may cause issues.

## Installation

### Conda environment

We recommend using a dedicated Conda environment for AbRaCD. This keeps Python and PyRosetta dependencies isolated from your other projects. Run the following commands from the repository root:

```bash
conda create -y -n abracd
conda activate abracd

conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta
python -m pip install mpi4py
bash install_path_conda.sh

conda deactivate
conda activate abracd
```

The PATH setup is scoped to the `abracd` environment. Deactivating and reactivating it applies the setup; AbRaCD commands are available while the environment is active.

### General Bash user without Conda

If you do not use Conda, you can add the repository's `scripts` and `src` command directories to your Bash PATH instead. Run this from the repository root:

```bash
python -m pip install mpi4py

bash install_path_bash.sh
```

Open a new terminal or run `source ~/.bashrc` to apply the change. This option makes the AbRaCD commands available in all Bash sessions for your user.
Install [PyRosetta](https://www.pyrosetta.org/downloads) separately in the Python environment you will use to run AbRaCD; it is required for structure cleaning and loop refinement.

### RCD

AbRaCD is integrated into [RCD](https://github.com/chaconlab/RCD.git). Follow the installation instructions in the RCD repository.

Here you have some secure steps to make everything work:

#### - System packages

Install the general build tools and dependencies for the GCC version of RCD if they are not already installed:

```bash
sudo apt update
sudo apt install -y build-essential cmake liblapack-dev
```

MPI is required to build and run an MPI version of RCD. Open MPI can be installed with:

```bash
sudo apt install -y openmpi-bin libopenmpi-dev
```

#### - Build RCD

Run the commands below from the root of the cloned RCD repository. We recommend building at least one MPI version. You can build any combination of the four variants; use a separate build directory for each one.

```bash
cmake -B build            && cmake --build build -j                  # single processor, GCC  -> bin/rcd_gcc
cmake -B build-mpi -DRCD_MPI=ON   && cmake --build build-mpi -j      # MPI                    -> bin/rcd_mpi_gcc
cmake -B build-intel -DRCD_INTEL=ON && cmake --build build-intel -j  # Intel icpx + MKL       -> bin/rcd_intel
cmake -B build-impi -DRCD_INTEL=ON -DRCD_MPI=ON                      # Intel + MPI            -> bin/rcd_mpi_intel
```

The GCC builds require LAPACK. MPI builds also require an MPI implementation and its development files. Intel builds require Intel oneAPI with `icpx` and MKL; load the oneAPI environment before building. The Intel + MPI build also requires a compatible MPI installation.

#### - Download the KORP energy file

The **KORP energy file** is not included in the RCD repository. From the root of the cloned RCD repository, download it with:

```bash
./download_korp.sh
```

The script saves `korp6Dv1.bin` to `RCD/src/libenergy/` by default. The Dunbrack file, `dunbrack.bin`, is included in that same directory.

#### - Make the RCD executables available on `PATH`

Add the `bin` directory of the cloned RCD repository to `PATH` so the executables can be called from any working directory. Replace `/absolute/path/to/RCD` with the path to your RCD clone:

```bash
export PATH="/absolute/path/to/RCD/bin:$PATH"
```

This applies to the current terminal. To make it persistent, add the line to `~/.bashrc` for general Bash use, or to a Conda activation script in the environment used for AbRaCD. Open a new terminal or reactivate the Conda environment after saving the change.

#### - Copy the RCD energy files to AbRaCD

AbRaCD expects the KORP and Dunbrack files in `data/Maps/`. Set `RCD_DIR` to the root of the cloned RCD repository and run these commands from the AbRaCD repository root:

```bash
RCD_DIR="/absolute/path/to/RCD"
cp "$RCD_DIR/src/libenergy/korp6Dv1.bin" data/Maps/
cp "$RCD_DIR/src/libenergy/dunbrack.bin" data/Maps/
```

## Quick start

Run commands from the directory where you want generated inputs, results, and logs. The main entry point is `abracd.pl`, available after running `install_path_conda.sh` or `install_path_bash.sh`. Required maps must be in `data/Maps/`, and an RCD executable must be available on `PATH`. See [AbRaCD_Manual.md](AbRaCD_Manual.md) for setup checks, input formats, options, and outputs.

```bash
# Show the command help
abracd.pl --help

# Chothia-numbered PDB: preprocess and model H3
abracd.pl antibody.pdb --CDRs H3

# Select several CDRs, or all six
abracd.pl antibody.pdb --CDRs H1 H3 L1
abracd.pl antibody.pdb --CDRs ALL

# Explicitly select source heavy/light chains using a purge table
abracd.pl purge.txt --purge --CDRs H3

# Run from an existing RCD input list
abracd.pl AbRaCD_input_H3.txt --raw --CDRs H3

# Design/modeling workflow from an RCD input list
abracd.pl AbRaCD_GAP_H3.txt --model --CDRs H3

# Configure MPI, sampling, repetitions, and an output/log prefix
abracd.pl antibody.pdb --CDRs H3 --NP 8 --NLOOPS 50000 --LOCOBEST 1000 --repRCD 2 --repROS 1 --OUT run_h3

# Run RCD without the Rosetta refinement stage
abracd.pl antibody.pdb --CDRs H3 --repROS 0 --OUT rcd_only_h3
```

The run log is `AbRaCD.log`, or `<OUT>.log` when `--OUT` is set. `--NP` enables MPI RCD execution; omit it for serial RCD.

## Fast commands for reproducing benchmark paper results

### AbRaCD Leakege Free Benchmark (ALFB)

The repository includes preprocessed ALFB benchmark inputs in `data/Benchmarks/ALFB/Processed/`, alongside the renumbered PDB structures. Since all six files follow the `AbRaCD_input_<CDR>.txt` naming convention expected by raw mode, `--CDRs ALL` runs the complete set in one invocation:

```bash
cd data/Benchmarks/ALFB/Processed
abracd.pl --raw --CDRs ALL --NP 8 --OUT ALFB
```

Runs all six CDRs with the default sampling settings and Rosetta refinement. Adjust `--NP` to match your available MPI processes. Results use the `ALFB` prefix; the log is `ALFB.log`.

### RosettaBenchmark

The Rosetta benchmark also has all six prepared `AbRaCD_input_<CDR>.txt` lists in its `Processed/` directory. Run the six CDRs together with `ALL`:

```bash
cd data/Benchmarks/RosettaBenchmark/Processed
abracd.pl --raw --CDRs ALL --NP 8 --OUT RBench
```

Outputs use the `RBench` prefix and the log is `RBench.log`.

## Repository contents

- `scripts/perl/abracd.pl`: end-to-end command-line workflow.
- `scripts/perl/AbRaCD_chothia_preprocces.pl`: Chothia preprocessing and RCD input generation.
- `scripts/perl/ab_engineering_prep.pl`: design structure preparation used by `--model`.
- `scripts/perl/refine_loops.pl` and `src/python/`: Rosetta refinement and structure utilities.
- `data/Maps/`: CDR, KORP, and Dunbrack maps used by RCD.
- `data/Benchmarks/`: benchmark structures and example RCD inputs.

## Troubleshooting

Run `abracd.pl --help` for the current syntax. Check the log for the command and stage that failed. Frequent setup issues include missing map files, missing RCD executables on `PATH`, or a Python environment without PyRosetta. See [AbRaCD_Manual.md](AbRaCD_Manual.md) for details.
