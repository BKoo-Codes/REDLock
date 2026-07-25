import os
import pandas as pd
from datetime import datetime
from .utils import calc_dg_nn, calculate_gc, tm_nn, rc
from .thermo import run_rnacofold

def print_split_report(res):
    report = f"""
[Target: {res['mirna_name']}]
  Optimal Split: {len(res['mirna_left'])} + {len(res['mirna_right'])} nt
  
  5' Arm: {res['mirna_left']} 
     (GC={res['gc_l']:.1f}%, Tm={res['tm_l']:.1f}°C) -> Padlock 3' end: {res['mirna_left_rc']}
     
  3' Arm: {res['mirna_right']} 
     (GC={res['gc_r']:.1f}%, Tm={res['tm_r']:.1f}°C) -> Padlock 5' end: {res['mirna_right_rc']}
  ------------------------------------------------------------"""
    print(report)

def _get_seq_analysis_details(title, seq):
    """Helper to generate a formatted string for a single sequence analysis."""
    if not seq or seq == "N/A":
        return ""
        
    gc_val = calculate_gc(seq)
    tm_val = tm_nn(seq)
    dg_val = calc_dg_nn(seq)
    
    report_title = f"===== [ {title} Analysis ] ====="
    return (
        f"{report_title}\n"
        f"  Sequence          : {seq}\n"
        f"  Length            : {len(seq)} nt\n"
        f"  GC Content        : {gc_val:.2f}%\n"
        f"  Melting Temp. (Tm): {f'{tm_val:.2f} °C' if tm_val is not None else 'N/A'}\n"
        f"  Folding ΔG        : {f'{dg_val:.2f} kcal/mol' if dg_val is not None else 'N/A'}\n"
        f"--------------------------------------------------\n"
    )

def _get_duplex_analysis_details(title, seq1, seq2, rnafold_path=None):
    """Helper to generate a more readable formatted string for a duplex analysis."""
    if not all([seq1, seq2, seq1 != "N/A", seq2 != "N/A"]):
        return ""

    _, duplex_dg = run_rnacofold(seq1, seq2)
    
    # --- Get names from title ---
    try:
        name1, name2 = title.split(' vs ')
    except ValueError:
        name1, name2 = "Seq1", "Seq2"

    # --- Properties for Seq1 ---
    gc1, tm1 = calculate_gc(seq1), tm_nn(seq1)
    # --- Properties for Seq2 ---
    gc2, tm2 = calculate_gc(seq2), tm_nn(seq2)

    report_title = f"===== [ Interaction: {title} ] ====="
    
    # Left-align names for clean layout
    name1_padded = f"{name1}:".ljust(18)
    name2_padded = f"{name2}:".ljust(18)

    return (
        f"{report_title}\n"
        f"  - {name1_padded} {seq1} (Len: {len(seq1)}, GC: {gc1:.1f}%, Tm: {tm1:.1f}°C)\n"
        f"  - {name2_padded} {seq2} (Len: {len(seq2)}, GC: {gc2:.1f}%, Tm: {tm2:.1f}°C)\n\n"
        f"  > Predicted Binding Energy (ΔG): {duplex_dg if duplex_dg is not None else 'N/A'} kcal/mol\n"
        f"--------------------------------------------------\n"
    )

def generate_detailed_report(row):
    """
    Generates a comprehensive, formatted string report for a single padlock candidate.
    """
    # Extract sequences, handling potential missing data
    padlock_full = row.get("Padlock Sequence (5'-P)", "N/A")
    padlock_seq = padlock_full.replace("[5'-P]-", "")
    spacer_seq = row.get("Spacer (Pad)", "N/A")
    cas12a_core_seq = row.get("Cas12a Core (Pad)", "N/A")
    
    mirna_seq = row.get("miRNA Sequence", "N/A")
    mirna_left_rc = row.get("miRNA Left Arm (RC)", "N/A")
    mirna_right_rc = row.get("miRNA Right Arm (RC)", "N/A")
    
    rca_product_seq = rc(padlock_seq)

    # Start building the report string
    report = ""

    # --- Single Sequence Analyses ---
    report += _get_seq_analysis_details("Full Padlock", padlock_seq)
    report += _get_seq_analysis_details("Spacer", spacer_seq)
    report += _get_seq_analysis_details("Cas12a Core", cas12a_core_seq)
    
    report += "\n" # Spacer

    # --- Duplex Analyses ---
    report += _get_duplex_analysis_details("Padlock vs Cas12a", padlock_seq, cas12a_core_seq)
    report += _get_duplex_analysis_details("Padlock vs miRNA left", padlock_seq, mirna_left_rc)
    report += _get_duplex_analysis_details("Padlock vs miRNA right", padlock_seq, mirna_right_rc)
    report += _get_duplex_analysis_details("RCA Product vs Cas12a", rca_product_seq, cas12a_core_seq)
    report += _get_duplex_analysis_details("RCA Product vs miRNA Sequence", rca_product_seq, mirna_seq)

    return report


def print_detailed_analysis(title, seq):
    """
    Comprehensive analysis of sequence information (GC, Tm, ΔG).
    """
    gc_val = calculate_gc(seq)
    tm_val = tm_nn(seq)
    dg_val = calc_dg_nn(seq)

    print(f"  [{title}]")
    print(f"    Seq    : {seq}")
    print(f"    Length : {len(seq)} nt")
    print(f"    GC%    : {gc_val:.2f}%")
    print(f"    Tm     : {f'{tm_val:.2f} °C' if tm_val is not None else 'N/A'}")
    print(f"    ΔG     : {f'{dg_val:.2f} kcal/mol' if dg_val is not None else 'N/A'}")
    print("-" * 50)

def export_candidate_results(data_list, filename_prefix="validated_padlocks", target_dir="results"):
    """
    Exports the list of validated candidates to an Excel file with a timestamp.
    
    Args:
        data_list (list): List of dictionaries or Series objects containing candidate data.
        filename_prefix (str): Prefix for the output filename.
        target_dir (str): Directory where the file will be saved.
    """
    # 1. Validation Check
    if not data_list:
        print("\n[WARNING] EXPORT FAILED: NO CANDIDATES FOUND IN THE LIST.")
        print("[WARNING] PLEASE CHECK IF THE VALIDATION LOOP COMPLETED SUCCESSFULLY.")
        return None

    # 2. Directory Setup
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
        print(f"[INFO] CREATED DIRECTORY: {target_dir}")

    # 3. Data Processing
    final_df = pd.DataFrame(data_list)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_filename = f"{filename_prefix}_{timestamp}.xlsx"
    file_path = os.path.join(target_dir, output_filename)

    # 4. Save to Excel
    try:
        final_df.to_excel(file_path, index=False)
        print(f"\n[SUCCESS] {len(final_df)} CANDIDATES EXPORTED TO EXCEL.")
        print(f"[SUCCESS] PATH: {file_path}")
        return final_df
    except Exception as e:
        print(f"\n[ERROR] FAILED TO SAVE EXCEL FILE: {e}")
        return None