import os
import pandas as pd
from datetime import datetime
from .utils import calc_dg_nn, calculate_gc, tm_nn

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