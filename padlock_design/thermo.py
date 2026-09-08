import glob
import os
import subprocess
import re
from shutil import which


class RNAfoldError(RuntimeError):
    """Raised when folding ΔG cannot be computed.

    Previously run_rnafold_dna() returned ("N/A", 0.0) on failure. Because
    ΔG = 0.0 is *above* every sane stability threshold, a missing RNAfold made
    every candidate silently PASS with a bogus penalty score. Failing loudly
    matches the documented behaviour in README ("the tool will stop rather than
    silently return ΔG = 0").
    """


class RNAfoldNotFoundError(RNAfoldError):
    """RNAfold executable could not be located."""


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

# Directories to sweep for dna_mathews*.par when none of the explicit Windows
# candidates exist. Needed because RNAfold exits with an error when -P points at
# a missing file — which, now that failures are no longer swallowed, would stop
# every non-Windows run.
_DNA_PAR_SEARCH_DIRS = [
    "/usr/share/ViennaRNA/misc",
    "/usr/local/share/ViennaRNA/misc",
    "/opt/homebrew/share/ViennaRNA/misc",
    "/opt/local/share/ViennaRNA/misc",
]


def _resolve_dna_par(executable=None):
    """Finds a ViennaRNA DNA parameter file (dna_mathews*.par).

    Search order: the VIENNA_DNA_PARAMS environment variable, the explicit
    (Windows) candidates above, directories relative to the RNAfold executable,
    then the standard Unix install locations. Returns "" when nothing is found.
    """
    env_path = os.environ.get("VIENNA_DNA_PARAMS", "")
    if env_path and os.path.isfile(env_path):
        return env_path

    for path in DNA_PAR_CANDIDATES:
        if os.path.isfile(path):
            return path

    search_dirs = list(_DNA_PAR_SEARCH_DIRS)
    if executable:
        exe_dir = os.path.dirname(os.path.abspath(executable))
        search_dirs[:0] = [
            exe_dir,
            os.path.join(exe_dir, "misc"),
            os.path.join(exe_dir, "..", "share", "ViennaRNA", "misc"),
            os.path.join(exe_dir, "..", "misc"),
        ]
    search_dirs.append(os.path.dirname(os.path.abspath(__file__)))

    for directory in search_dirs:
        for match in sorted(glob.glob(os.path.join(directory, "dna_mathews*.par"))):
            if os.path.isfile(match):
                return os.path.normpath(match)
    return ""


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

    Raises:
        RNAfoldNotFoundError: the executable could not be located.
        RNAfoldError: RNAfold ran but produced unusable output.
    """
    executable = rnafold_executable if rnafold_executable else which("RNAfold")
    if not executable:
        raise RNAfoldNotFoundError(
            "RNAfold executable not found. Install ViennaRNA and put RNAfold on "
            "your PATH, or pass an explicit path (--rnafold_path / the RNAfold "
            "Path field in the GUI). Folding ΔG cannot be computed without it."
        )
    if rnafold_executable and not (os.path.isfile(executable) or which(executable)):
        raise RNAfoldNotFoundError(f"RNAfold executable not found at: {executable}")

    # Re-resolve using the executable as a hint (handles Linux/macOS installs).
    resolved_dna_par_path = (RESOLVED_DNA_PAR if os.path.isfile(RESOLVED_DNA_PAR)
                             else _resolve_dna_par(executable))
    if not resolved_dna_par_path:
        raise RNAfoldError(
            "ViennaRNA DNA parameter file (dna_mathews*.par) not found. DNA "
            "parameters are required for padlock folding ΔG; pointing RNAfold at "
            "a missing file makes it exit with an error. Set the "
            "VIENNA_DNA_PARAMS environment variable to the file's full path."
        )

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
    except OSError as e:
        raise RNAfoldError(f"RNAfold execution failed ({executable}): {e}") from e

    if process.returncode != 0:
        raise RNAfoldError(
            f"RNAfold exited with code {process.returncode}: {(stderr or '').strip()}"
        )

    lines = stdout.strip().splitlines()
    if len(lines) < 2:
        raise RNAfoldError(
            f"Unexpected RNAfold output for sequence of length {len(seq)}: {stdout!r}"
        )

    structure = lines[1].split()[0]
    # Energy is the last parenthesised number on the line, e.g. "....(( ))  ( -3.20)"
    match = re.search(r"\(\s*([-+]?\d+(?:\.\d+)?)\s*\)\s*$", lines[1])
    if not match:
        raise RNAfoldError(f"Could not parse MFE from RNAfold output line: {lines[1]!r}")

    return structure, float(match.group(1))

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
        par = RESOLVED_DNA_PAR if os.path.isfile(RESOLVED_DNA_PAR) else _resolve_dna_par(exe)
        if par:
            cmd.extend(["--noconv", "-P", par])

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