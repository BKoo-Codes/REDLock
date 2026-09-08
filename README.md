# Padlock Probe Design Tool

REDLock is an automated design pipeline for creating high-specificity padlock probes for miRNA detection using RCA (Rolling Circle Amplification) and CRISPR-Cas12a systems.

## Note

Personal working copy — kept for development and record-keeping.
------------------------------------------------------------------------------------

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/BKoo-Codes/REDLock.git
   cd REDLock
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

   *Note: You may need to install NCBI BLAST+ locally if using local BLAST features, although this tool currently uses the NCBI web API via Biopython.*

## Features

- **Automated Design**: Generates optimized spacer and Cas12a core sequences.
- **Strict Filtering**: Checks for homopolymers, hairpins, and cross-hybridization risks.
- **Penalty-Based Ranking**: Scores candidates on GC content, folding ΔG, and arm / Cas12a accessibility.
- **Stability Analysis**: Calculates ΔG and predicts secondary structures using ViennaRNA (DNA parameters).
- **NCBI Validation**: BLAST validation against the human transcriptome, including contiguous-match and 3'-terminal-match checks.

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd REDLock
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Install [ViennaRNA](https://www.tbi.univie.ac.at/RNA/). `RNAfold` must be available,
   either on your `PATH`, next to the application, or at a standard install location.
   Without it, folding ΔG cannot be computed and the tool will stop rather than
   silently return ΔG = 0.

   Folding uses ViennaRNA's DNA parameters (`dna_mathews*.par`), which ship with
   ViennaRNA. The file is looked up automatically next to the `RNAfold`
   executable and at the usual install locations; if yours lives elsewhere, set
   the `VIENNA_DNA_PARAMS` environment variable to its full path.

## Pool cache

Spacer / Cas12a pools are cached to a `.pkl` file. The cache records the
parameters it was built from (lengths, GC ranges, sample count, seed) and is
reused only when all of them match the current run, so changing a setting always
regenerates the pools instead of silently reusing stale ones.

## Usage

### Option 1: Graphical User Interface (GUI) - Recommended

Launch the interactive web-based dashboard:

```bash
streamlit run app.py
```

This will open the design tool in your default web browser. You can configure parameters, run the design pipeline, and download results directly from the interface.

### Option 2: Command Line Interface (CLI)

Run the tool via the terminal for batch processing or automation:

```bash
python main.py --target_name "hsa-miR-141-3p" \
               --target_seq "UAACACUGUCUGGUAAAGAUGG" \
               --num_candidates 2000 \
               --skip_blast \
               --output "results/my_design.xlsx"
```

### Option 3: Create Standalone Executable (.exe)

You can package this tool into a single executable file that runs without installing Python.

1.  **Install PyInstaller:**
    ```bash
    pip install pyinstaller
    ```

2.  **Build the Executable:**
    ```bash
    pyinstaller --onefile ^
      --add-data "app.py;." ^
      --add-data "padlock_design;padlock_design" ^
      --add-binary "C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe;." ^
      --name REDLock_Tool --clean run_gui.py
    ```
    *(On Linux/macOS, replace `;` with `:` in the `--add-data` / `--add-binary` arguments)*

    Bundling `RNAfold` is required for penalty scoring to work on machines
    without a ViennaRNA installation.

3.  **Run:**
    - Go to the `dist` folder and double-click `REDLock_Tool.exe`.

### Arguments (CLI)

| Argument | Default | Description |
| --- | --- | --- |
| `--target_name` | `hsa-miR-21-5p` | Name of the target miRNA |
| `--target_seq` | `UAGCUUAUCAGACUGAUGUUGA` | Sequence of the target miRNA |
| `--num_candidates` | `10` | Number of candidates to assemble before analysis |
| `--output` | `results/padlock_results.xlsx` | Path to save the Excel result file |
| `--skip_blast` | off | Skip BLAST validation (much faster) |
| `--dg_threshold` | `-5.0` | Candidates with ΔG at or below this are rejected as too stable |
| `--rnafold_path` | auto-detect | Path to the ViennaRNA `RNAfold` executable |
| `--spacer_length` | `10` | Spacer length to prioritize |
| `--cas12a_length` | `20` | Cas12a core length to prioritize |
| `--pool_seed` | `15` | Seeds pool generation |
| `--seed_offset` | `0` | Assembly seed start (distinct from the pool seed) |
| `--pool_samples` | `500000` | Random samples per pool length |
| `--blast_count` | `5` | Number of candidates to keep after BLAST validation |
| `--pool_file` | auto | Pool cache file. Defaults to a name derived from the pool parameters |


## Website

Visit our [Project Website](https://bkoo-codes.github.io/REDLock/) for more details.



## License

This project is licensed under the MIT License.
