# PadDx_miR: Padlock Probe Design Tool

PadDx_miR is an automated design pipeline for creating high-specificity padlock probes for miRNA detection using RCA (Rolling Circle Amplification) and CRISPR-Cas12a systems.

## Features

- **Automated Design**: Generates optimized spacer and Cas12a core sequences.
- **Strict Filtering**: Checks for homopolymers, hairpins, and cross-hybridization risks.
- **Stability Analysis**: Calculates ΔG and predicts secondary structures using RNAfold logic.
- **NCBI Validation**: Integrated BLAST validation to ensure specificity against the human transcriptome.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/BKoo-Codes/PadDx_miR.git
   cd PadDx_miR
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

   *Note: You may need to install NCBI BLAST+ locally if using local BLAST features, although this tool currently uses the NCBI web API via Biopython.*

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
python main.py --target_name "hsa-miR-21-5p" \
               --target_seq "UAGCUUAUCAGACUGAUGUUGA" \
               --num_candidates 10 \
               --output "results/my_design.xlsx"
```

### Option 3: Create Standalone Executable (.exe)

You can package this tool into a single executable file that runs without installing Python (ideal for sharing with non-technical users).

1.  **Install PyInstaller:**
    ```bash
    pip install pyinstaller
    ```

2.  **Build the Executable:**
    Run the following command in your terminal (Windows PowerShell or Command Prompt):
    ```bash
    pyinstaller --onefile --add-data "app.py;." --add-data "padlock_design;padlock_design" --name PadDx_miR_Tool --clean run_gui.py
    ```
    *(Note: On Linux/macOS, replace `;` with `:` in the `--add-data` arguments)*

3.  **Run:**
    - Go to the `dist` folder.
    - You will find `PadDx_miR_Tool.exe`.
    - Double-click to run the application instantly!

### Arguments (CLI)

- `--target_name`: Name of the target miRNA (default: "hsa-miR-21-5p")
- `--target_seq`: Sequence of the target miRNA (default: "UAGCUUAUCAGACUGAUGUUGA")
- `--num_candidates`: Number of candidates to generate (default: 10)
- `--output`: Path to save the Excel result file (default: "results/padlock_results.xlsx")
- `--skip_blast`: Skip the NCBI BLAST validation step (faster execution)
- `--rnafold_path`: Path to the ViennaRNA `RNAfold` executable (e.g., `C:\\Program Files (x86)\\ViennaRNA Package\\RNAfold.exe` on Windows, or `/usr/local/bin/RNAfold` on Linux/macOS if not in PATH). If not provided, the tool will try to find `RNAfold` in your system's PATH.


## Website

Visit our [Project Website](https://bkoo-codes.github.io/PadDx_miR/) for more details.

## License

This project is licensed under the MIT License.
