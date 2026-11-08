# AbRaCD command-line manual

This manual describes the command-line interface implemented by `scripts/perl/abracd.pl` and the input files it consumes. Run commands from the working directory where you want the generated files, result folders, and log. The launcher can be checked with `abracd.pl --help`.

## 1. Requirements and setup checks

Use the installation steps in [Readme.md](Readme.md). Before a run, confirm that:

- `abracd.pl` is on `PATH` and Perl can execute it.
- Python can import `pyrosetta`. 
- MPI runs also require `mpi4py` and `mpirun`.
- `data/Maps/CDRS.bin`, `data/Maps/korp6Dv1.bin`, and `data/Maps/dunbrack.bin` exist.
- An RCD executable is on `PATH`. Serial selection looks for `rcd_gcc` or `rcd`; MPI runs look for `rcd_mpi_gcc`, `rcd_mpi_intel`, or `rcd_mpi` as well as `mpirun`. If only an MPI executable is available, the launcher can fall back to it.
- Rosetta refinement helpers are available through the project installation. `--resROS 0` skips Rosetta refinement.

The launcher checks several dependencies and input files before starting. Its log records the invocation and executed commands.

## 2. Main command

```text
abracd.pl [INPUT] [--raw | --model] [--purge] --CDRs CDR [CDR ... | ALL] [--OUT NAME] [--NP N] [--NLOOPS N] [--LOCOBEST N] [--repRCD N] [--repROS N]
```

`INPUT` is optional except where a mode requires it. `--CDRs` is required. Options are case-sensitive. If used, `--raw` or `--model` must immediately follow `INPUT` (or come first when there is no input). The selected CDRs run in the order given. `ALL` expands to `H1 H2 H3 L1 L2 L3` and must be used alone.

In the default mode, the command runs three stages:

1. **Chothia preprocessing:** prepares and renumbers the antibody structures and generates the RCD input lists for the selected CDRs. Use `--purge` in this stage when source chains must be explicitly mapped to heavy (H) and light (L); this option is only for Chothia-formatted structures.
2. **RCD sampling and selection:** generates loop decoys and selects candidates using the RCD method.
3. **Rosetta refinement:** refines the RCD results with PyRosetta.

See [Input modes](#input-modes) for the input requirements for each mode and [File formats](#3-file-formats) for the structure and RCD list formats.

### Input modes

| Mode | Input and behavior |
|---|---|
| Default preprocessing | `INPUT` may be one `.pdb` structure or a `.txt` list of PDB paths. With no input, all original `*.pdb` files in the current directory are processed. Chothia preprocessing runs before RCD for the requested CDRs. |
| `--purge` | Provide the purge table as `INPUT`; it selects structures and source H/L chain IDs. Only valid with default Chothia preprocessing. The PDB structures must use Chothia numbering. See [Purge table](#purge-table) for the required format. |
| `--raw` | Skips Chothia preprocessing. For one CDR, provide that CDR's RCD-format input list as `INPUT`. For multiple CDRs or `ALL`, omit `INPUT`; the launcher reads `AbRaCD_input_<CDR>.txt` files from the current directory. Since this mode does not prepare or validate the structures, ensure each PDB and RCD list is correct before running. As a practical guideline, keep only the two antibody chains, named `H` and `L`, and verify the loop start and end residue indices; incorrect chain or loop definitions can produce invalid results or cause the run to stall. See [RCD input list](#rcd-input-list). |
| `--model` | Design/modeling mode without native-structure RMSD. Provide one RCD-format `.txt` input and exactly one CDR. The launcher prepares design PDBs, runs RCD and refinement, and removes native-dependent score columns. Cannot be combined with `--raw` or `--purge`. See [Engineering preparation](#engineering-preparation) for the design preparation step and input requirements. |

### Options

| Option | Value | Required or default | Description |
|---|---|---|---|
| `--CDRs` | One or more of `H1 H2 H3 L1 L2 L3`, or `ALL` | **Required** | Selects which CDRs to process. CDRs run in the order listed. Duplicates are not allowed. `ALL` expands to all six CDRs and must be used by itself. |
| `--purge` | None | Optional; off by default | Treats `INPUT` as a purge table. The table selects structures and maps their source heavy and light chain IDs to `H` and `L`. Only applies to default Chothia preprocessing; it cannot be combined with `--raw` or `--model`. |
| `--raw` | None | Optional; off by default | Skips Chothia preprocessing and uses prepared RCD input lists. With one CDR, `INPUT` must be that CDR's RCD list. With multiple CDRs or `ALL`, omit `INPUT` and put the required `AbRaCD_input_<CDR>.txt` files in the current directory. |
| `--model` | None | Optional; off by default | Runs the design/modeling workflow. Requires one RCD-format `.txt` input and exactly one CDR. Prepares the design structure, runs RCD and Rosetta refinement in modeling mode, and removes native-dependent RMSD data. Cannot be combined with `--raw` or `--purge`. |
| `--OUT` | `NAME` | Optional; no prefix by default | Adds a prefix to generated RCD output names and uses `NAME.log` for the run log. `NAME` may contain letters, digits, periods, underscores, and hyphens. |
| `--NP` | Positive integer | Optional; serial RCD by default | Requests MPI execution of RCD with the specified number of processes. Requires an MPI RCD executable and `mpirun`. Without `--NP`, the launcher prefers serial RCD; if only an MPI executable is available, it falls back to MPI with one process. |
| `--NLOOPS` | Positive integer | Optional; `50000` for H3, `20000` for other CDRs | Overrides the number of loop samples requested from RCD for each selected CDR. |
| `--LOCOBEST` | Positive integer | Optional; `1000` for H3, `2000` for other CDRs | Overrides the RCD `loco_best` selection value for each selected CDR. |
| `--repRCD` | Positive integer | Optional; `1` | Sets the number of RCD sampling runs per selected CDR. |
| `--repROS` | Non-negative integer | Optional; `1` | Sets the number of Rosetta refinement runs for each RCD run. Use `0` to skip the refinement stage. |

Each option may be specified at most once. Integer parameters must be positive, except `--repROS`, which accepts zero. `--purge`, `--raw`, and `--model` are flags and take no value. The Rosetta repetition option is named `--repROS`.

### Examples

```bash
# 1. Run one Chothia-numbered antibody structure for one CDR.
abracd.pl antibody.pdb --CDRs H3

# 2. Process a text list of structures for several selected CDRs.
abracd.pl antibodies.txt --CDRs H1 H3 L1

# 3. Process all six CDRs for original PDB files in the current directory.
abracd.pl --CDRs ALL

# 4. Use a purge table to map the source chains to heavy (H) and light (L).
#    The structures must be Chothia-numbered. See the purge table format below.
abracd.pl purge.txt --purge --CDRs H3

# 5. Request MPI RCD execution with eight processes.
abracd.pl antibody.pdb --CDRs H3 --NP 8

# 6. Run RCD on a prepared input list without Chothia preprocessing.
#    Check the PDB chains and loop start/end indices before using --raw.
abracd.pl AbRaCD_input_H3.txt --raw --CDRs H3

# 7. Override both RCD sampling and candidate-selection counts.
abracd.pl antibody.pdb --CDRs H3 --NLOOPS 30000 --LOCOBEST 1500

# 8. Run the design/modeling workflow for one CDR.
#    See Engineering preparation for the design input and preparation details.
abracd.pl AbRaCD_GAP_H3.txt --model --CDRs H3

# 9. Run three RCD sampling repetitions and two Rosetta refinements per RCD run.
abracd.pl antibody.pdb --CDRs H3 --repRCD 3 --repROS 2

# 10. Run RCD but skip the Rosetta refinement stage.
abracd.pl antibody.pdb --CDRs H3 --repROS 0
```

Invalid combinations include `--raw` and `--model` together, `--purge` with either mode, a positional input for `--raw` with multiple CDRs, and multiple CDRs with `--model`. The mode flag must immediately follow the optional input, for example `input.txt --model --CDRs H3`.

## 3. File formats

### PDB list

A `.txt` structure list contains one PDB path per line. Relative paths are resolved from the directory containing the list. Blank lines are ignored. In default preprocessing mode, use either one PDB directly, or this list.

```text
structures/antibody_a.pdb
structures/antibody_b.pdb
```

When no input is provided, the preprocessor scans the current directory for original `.pdb` files (excluding files ending in `_renum.pdb`).

### Purge table

A purge table contains one whitespace-separated row per structure:

```text
PDB_path  H_chain_ID  L_chain_ID
```

Relative PDB paths are resolved from the purge file's directory. Use `-` for a missing chain, such as a nanobody. Blank lines and lines beginning with `#` are ignored. The structure name may be written with or without `.pdb` in the mapping lookup. The table defines which source chains are retained and mapped to H/L.

```text
# structure heavy-chain light-chain
complex_a.pdb A B
complex_b B Z
nanobody.pdb H -
```

Run it through the main workflow as `abracd.pl purge.txt --purge --CDRs H3`, or invoke preprocessing directly as described below.

### RCD input list

Each non-comment row has at least seven whitespace-separated fields:

```text
PDB_file  start  end  chain  sequence  length  CDR
```

`start` and `end` are inclusive residue indices; `end` must be greater than or equal to `start`. The sequence length must equal `length`. Chain is a one-character chain ID. Relative PDB paths are resolved from the RCD list directory (the design-preparation helper also accepts a matching path in the current directory).

```text
antibody_renum.pdb  97  109  H  ARWGGDGFYAMDY  13  H3
```

Comments and blank lines are ignored by the design preparation helper. Standard generated files are named `AbRaCD_input_<CDR>.txt`. For `--model`, the input is expected to describe the loop/design interval in this same RCD row format; the helper writes `AbRaCD_design_<CDR>.txt` and `<source>_design.pdb` beside the input structures/list.

## 4. Chothia preprocessing with `AbRaCD_chothia_preprocces.pl`

This helper prepares Chothia-numbered antibody structures for RCD. It selects and renumbers the antibody chains, identifies the CDR intervals, and writes the PDB structures and RCD input lists beside the source structures or input list. The main `abracd.pl` command runs this helper automatically in the default mode; it can also be run directly for preparation only.

```text
AbRaCD_chothia_preprocces.pl [--purge PURGE.txt] [PDB | PDB_LIST.txt]
```

Examples:

```bash
AbRaCD_chothia_preprocces.pl antibody.pdb
AbRaCD_chothia_preprocces.pl structures.txt
AbRaCD_chothia_preprocces.pl --purge purge.txt
AbRaCD_chothia_preprocces.pl --help
```

The input may be one PDB, a text list of PDB paths, or omitted to process original PDB files in the current directory. With `--purge`, provide the purge table described in [Purge table](#purge-table); it maps each source heavy and light chain ID to `H` and `L`. Without a purge table, accepted chain layouts are `H` and `L`, `H`/`L` plus one antigen chain, or one `H` or `L` antibody chain plus one antigen chain. Other unselected chains are removed. Residues beyond original Chothia H113 or L109 are removed before renumbering.

For each accepted structure, the helper writes a renumbered PDB named `<source>_renum.pdb` and generates `AbRaCD_input_<CDR>.txt` lists for the CDRs it identifies. If a critical gap is detected in a CDR, it can also write `AbRaCD_GAP_<CDR>.txt` for the design workflow described in [Section 5](#5-design-preparation-with-ab_engineering_preppl). The generated files are placed beside the input PDB or in the directory containing the PDB list.

## 5. Design preparation with `ab_engineering_prep.pl`

### What the preparation script does

`ab_engineering_prep.pl` prepares an antibody structure and its RCD input for modeling a loop with a user-specified sequence. It reads an RCD-format text file and, for each row, removes the residues in the specified inclusive start-to-end interval from the specified chain. It then shifts the numbering of downstream residues in that chain so the replacement loop can have the requested length. Other chains and non-coordinate PDB records are retained. The script does not build the new loop: it creates the structure and RCD list that the subsequent RCD and Rosetta stages use to generate and refine loop models.

The input list must contain the fields described in [RCD input list](#rcd-input-list): PDB path, start, end, chain, target sequence, target sequence length, and CDR. The sequence length must match the number of residues in the target sequence. Start and end identify the residues to remove from the input structure; ensure that the PDB path, residue numbering, chain, and interval all refer to the same structure.

```text
ab_engineering_prep.pl RCD_INPUT.txt
```

For example, to run only the preparation step:

```bash
ab_engineering_prep.pl design_input.txt
```

The script writes `<source>_design.pdb` and `AbRaCD_design_<CDR>.txt` beside the source structures and input list. If the source PDB is named `<source>_renum.pdb`, the output structure is named `<source>_design.pdb`. Running this helper alone stops after preparing these files; it does not run RCD or Rosetta.

The main command's `--model` mode calls this helper automatically, then passes the generated `AbRaCD_design_<CDR>.txt` list to RCD and Rosetta. For a complete modeling run, pass the original, user-prepared input list to `abracd.pl` with `--model` and the matching single CDR. Do not pass the generated `AbRaCD_design_<CDR>.txt` list: `--model` would treat it as a fresh preparation input and remove the loop interval again.

### Recommended workflow for Chothia-numbered structures

1. Run Chothia preprocessing on the source PDB or PDB list, either directly or through the default `abracd.pl` workflow. For preparation only:

   ```bash
   AbRaCD_chothia_preprocces.pl antibody.pdb
   ```

   If a structural gap is detected in a CDR, the preprocessor reports it and writes `AbRaCD_GAP_<CDR>.txt`. When running `abracd.pl` for a CDR with a detected critical gap, preprocessing reports the gap and stops before RCD so the design input can be reviewed.

2. Edit the relevant `AbRaCD_GAP_<CDR>.txt` row. The generated sequence and length describe the residues observed in the CDR, which may be incomplete. Replace them with the desired full loop sequence and its length. Check that the referenced `<source>_renum.pdb`, chain, and inclusive start/end interval identify the existing residues that should be replaced. Keep the sequence and length consistent.

3. Pass the edited GAP file to `abracd.pl --model`:

   ```bash
   abracd.pl AbRaCD_GAP_H3.txt --model --CDRs H3
   ```

   `abracd.pl` runs `ab_engineering_prep.pl` automatically, renumbers residues after the replaced interval to accommodate the target length, and uses the generated design input for RCD and Rosetta. There is no need to run the helper separately for this complete workflow.

### Recommended workflow for raw structures

Raw structures have not been prepared or renumbered by the Chothia preprocessor. Create an RCD-format input file yourself and ensure that its PDB path, chain ID, and residue numbering match the raw PDB. The input row must provide the replacement interval and target loop information:

```text
PDB_file  start  end  chain  target_sequence  target_length  CDR
```

The start and end values are inclusive and refer to residue numbers in the raw PDB. The script uses the integer residue numbers and chain ID to remove existing residues in this interval; make sure the interval covers the residues to replace. It does not distinguish PDB insertion codes, so use an interval whose integer residue numbers identify the intended residues unambiguously. Set `target_sequence` to the desired loop sequence and `target_length` to its exact length. This is more than adjusting the start and end fields: the target sequence and length must also be supplied. The helper does not infer a missing loop sequence from a raw structure.

Run the complete raw-structure design workflow by passing this original, user-prepared list to `abracd.pl --model`:

```bash
abracd.pl design_input.txt --model --CDRs H3
```

As with Chothia inputs, `abracd.pl --model` invokes the preparation helper and then runs RCD and Rosetta using the generated `AbRaCD_design_<CDR>.txt` file. You may run `ab_engineering_prep.pl design_input.txt` on its own to prepare and inspect the design structure and list, but a subsequent `abracd.pl --model` run should still receive the original `design_input.txt`; it will repeat the preparation. Do not pass the generated design list as the `--model` input.

## 6. Outputs and logs

The main command creates outputs for each selected CDR and RCD repetition. By default, the RCD result directory is named `<CDR>_rep<N>` (for example, `H3_rep1`). When `--OUT NAME` is used, the prefix is added (`<NAME>_<CDR>_rep<N>`). These directories and the run log are written relative to the directory from which `abracd.pl` is launched.

### Preprocessing outputs

Chothia preprocessing writes the renumbered and cleaned `<source>_renum.pdb` structures beside the input structures, together with `AbRaCD_input_<CDR>.txt` files containing the detected loop intervals and sequences. For a critical gap, it may also write `AbRaCD_GAP_<CDR>.txt`, which can be used as the input to the design workflow in [Section 5](#5-design-preparation-with-ab_engineering_preppl).

### RCD outputs

Each RCD run creates a result directory containing the generated loop decoys and `results.txt`. The results file includes RCD run information and, when a native crystallographic structure is available, the benchmark summary used to assess KORP filtering performance. H3 inputs containing loops of seven residues or fewer are split internally; those runs may also produce `results_shorts.txt` for the short-loop subset.

### Rosetta refinement outputs for native structures

When refinement runs with a known native structure, it writes refined loop structures and score tables inside the corresponding RCD result directory. The per-case `*_rosettaS.txt` table is sorted by Rosetta energy, with the lowest-energy candidates first. Summary tables named `rosetta1...txt`, `rosetta2...txt`, `rosetta5...txt`, `rosetta10...txt`, `rosetta20...txt`, `rosetta50...txt`, and `rosetta100...txt` report RMSD results for the best candidate among the first 1, 2, 5, 10, 20, 50, or 100 candidates in that energy ranking. Thus, for example, the “best RMSD of the top 5” is the smallest RMSD to the native structure among the five lowest-energy refined candidates; it does not mean the five candidates with the lowest RMSD were ranked first. Smaller RMSD indicates a closer match to the crystallographic loop.

### Design/modeling outputs

In design/modeling mode there is no native reference structure, so native-dependent RMSD columns and summaries are not available. The refined score tables remain ordered by energy; use the ranked `*_rosettaS.txt` results to compare the generated models and choose candidates. The corresponding refined structures are written as PDB files in the same RCD result directory.

### Run log

The launcher writes a detailed log in the launch directory: `AbRaCD.log` by default, or `<NAME>.log` when `--OUT NAME` is supplied. It records the command, selected CDRs and settings, stage commands, output directories, and errors. Use it to identify which stage ran and to diagnose a failed run.
