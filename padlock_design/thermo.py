import os
import subprocess
import re
from shutil import which

# ===== ViennaRNA 설정 =====
VIENNA_PATH = r"C:\Program Files (x86)\ViennaRNA Package"
RNACOFOLD_PATH = os.path.join(VIENNA_PATH, "RNAcofold.exe")
RNAFOLD_PATH = os.path.join(VIENNA_PATH, "RNAfold.exe")
DNA_PARAMETERS = os.path.join(VIENNA_PATH, "misc", "dna_mathews.par")
VIENNA_TEMP_C = 37.0    # φ29 RCA Defaults
VIENNA_DANGLES = 2

DNA_PAR_CANDIDATES = [
    os.path.join(VIENNA_PATH, "misc", "dna_mathews.par"),
    os.path.join(VIENNA_PATH, "share", "ViennaRNA", "misc", "dna_mathews.par"),
    r"C:\Program Files\ViennaRNA Package\misc\dna_mathews.par",
    r"C:\Program Files (x86)\ViennaRNA Package\Misc\dna_mathews2004.par",
]

def _resolve_dna_par():
    """Finds the first valid DNA parameter file path from candidates."""
    for path in DNA_PAR_CANDIDATES:
        if os.path.isfile(path):
            return path
    return "dna_mathews.par"

RESOLVED_DNA_PAR = _resolve_dna_par()

def run_rnafold_dna(seq, temp=37.0, rnafold_executable=None):
    """
    Predicts DNA secondary structure and MFE using RNAfold.

    Args:
        seq (str): DNA sequence.
        temp (float): Simulation temperature (°C). Defaults to 37.0.
        rnafold_executable (str, optional): Path to the RNAfold executable. 
                                            If None, attempts to find it in PATH.

    Returns:
        tuple[str, float]: (Dot-bracket structure, MFE in kcal/mol).
                           Returns ("N/A", 0.0) if an error occurs or executable not found.
    """
    executable = rnafold_executable if rnafold_executable else which("RNAfold")
    if not executable:
        print("[ERROR] RNAfold executable not found. Please provide path or ensure it's in PATH.")
        return "N/A", 0.0

    # Fallback for DNA_PARAMETERS if RESOLVED_DNA_PAR is not a valid file path
    resolved_dna_par_path = RESOLVED_DNA_PAR if os.path.isfile(RESOLVED_DNA_PAR) else "dna_mathews.par"

    cmd = [
            executable, 
            "--noconv", 
            "-P", resolved_dna_par_path, 
            f"--temp={temp}", 
            f"-d{VIENNA_DANGLES}",
            "--noPS"
        ]
    
    try:
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=seq)
        
        lines = stdout.strip().splitlines()
        if len(lines) >= 2:
            structure = lines[1].split()[0]
            # Extract energy, robustly handle cases where split might not yield expected parts
            energy_str = lines[1].split('(')[-1].strip(')')
            try:
                energy = float(energy_str)
            except ValueError:
                energy = 0.0 # Default to 0.0 if energy parsing fails
            return structure, energy
            
    except Exception as e:
        print(f"[ERROR] RNAfold execution failed: {e}")
        return "N/A", 0.0

    return "N/A", 0.0

def run_rnacofold(seq1, seq2, temperature=VIENNA_TEMP_C, dangles=VIENNA_DANGLES, 
                  use_dna_params=True, return_structure=True):
    """
    Calculates dimer delta G using RNAcofold.

    Args:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.
        temperature (float): Simulation temperature.
        dangles (int): Dangle model (0, 1, 2).
        use_dna_params (bool): Use DNA parameters if True.
        return_structure (bool): Returns (structure, dG) if True.

    Returns:
        tuple | float | None: (str, float) or float depending on return_structure.
    """
    # 1. Resolve executable path
    exe = RNACOFOLD_PATH if os.path.isfile(RNACOFOLD_PATH) else which("RNAcofold")
    if not exe:
        return (None, None) if return_structure else None

    # 2. Build command
    cmd = [exe, "-T", str(temperature), f"-d{dangles}", "--noPS"]
    if use_dna_params:
        cmd.extend(["--noconv", "-P", RESOLVED_DNA_PAR])

    # 3. Execute subprocess
    input_str = f"{seq1}&{seq2}\n"
    try:
        proc = subprocess.run(cmd, input=input_str, capture_output=True, text=True, check=True)
        out = proc.stdout.strip().splitlines()
        
        if len(out) > 1:
            line = out[1]
            # Extract energy value inside parentheses using regex
            match = re.search(r"\(\s*([-+]?\d*\.\d+|\d+)\s*\)", line)
            energy = float(match.group(1)) if match else 0.0
            structure = line.split()[0] if " " in line else ""
            
            return (structure, energy) if return_structure else energy
            
    except (subprocess.CalledProcessError, Exception):
        return (None, None) if return_structure else None

    return (None, None) if return_structure else None