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

Run the tool via the command line interface (CLI):

```bash
python main.py --target_name "hsa-miR-21-5p" \
               --target_seq "UAGCUUAUCAGACUGAUGUUGA" \
               --num_candidates 10 \
               --output "results/my_design.xlsx"
```

### Arguments

- `--target_name`: Name of the target miRNA (default: "hsa-miR-21-5p")
- `--target_seq`: Sequence of the target miRNA (default: "UAGCUUAUCAGACUGAUGUUGA")
- `--num_candidates`: Number of candidates to generate (default: 10)
- `--output`: Path to save the Excel result file (default: "results/padlock_results.xlsx")
- `--skip_blast`: Skip the NCBI BLAST validation step (faster execution)
- `--dg_threshold`: Delta G threshold for stability filtering (default: -3.5)

## Website

Visit our [Project Website](https://bkoo-codes.github.io/PadDx_miR/) for more details.

## License

This project is licensed under the MIT License.
