import argparse
import pandas as pd
import os
import sys
from padlock_design.pipeline import (
    generate_pool, 
    assemble_and_filter_padlocks, 
    analyze_and_sort_padlocks, 
    run_blast_validation,
    export_candidate_results
)

def main():
    parser = argparse.ArgumentParser(description="PadDx_miR: Padlock Probe Design Tool")
    
    parser.add_argument("--target_name", type=str, default="hsa-miR-21-5p", help="Target miRNA name")
    parser.add_argument("--target_seq", type=str, default="UAGCUUAUCAGACUGAUGUUGA", help="Target miRNA sequence")
    parser.add_argument("--num_candidates", type=int, default=10, help="Number of candidates to generate")
    parser.add_argument("--output", type=str, default="results/padlock_results.xlsx", help="Output Excel file path")
    parser.add_argument("--skip_blast", action="store_true", help="Skip NCBI BLAST validation (faster)")
    parser.add_argument("--dg_threshold", type=float, default=-3.5, help="Delta G threshold for stability filtering")
    
    args = parser.parse_args()
    
    # Create results directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"\n{'='*50}")
    print(f"PadDx_miR Padlock Design Tool")
    print(f"{'='*50}")
    print(f"Target: {args.target_name}")
    print(f"Sequence: {args.target_seq}")
    print(f"Candidates to generate: {args.num_candidates}")
    print(f"Output: {args.output}")
    print(f"{'='*50}\n")

    # 1. Generate/Load Pools
    print(">>> Step 1: Initializing Sequence Pools...")
    spacer_length = [16, 15, 17]
    cas12a_length = [20, 21, 22]
    spacer_pools, cas12a_pools = generate_pool(spacer_length, cas12a_length)

    # 2. Assemble and Filter
    print("\n>>> Step 2: Assembling and Filtering Candidates...")
    candidates = assemble_and_filter_padlocks(
        args.target_name, 
        args.target_seq, 
        spacer_pools, 
        cas12a_pools, 
        num_to_generate=args.num_candidates
    )

    if not candidates:
        print("[ERROR] No candidates generated. Check parameters.")
        sys.exit(1)

    df = pd.DataFrame(candidates)
    
    # 3. Analyze and Sort by Stability
    print("\n>>> Step 3: Analyzing Stability (RNAfold)...")
    sorted_df = analyze_and_sort_padlocks(df, dg_threshold=args.dg_threshold)
    
    if sorted_df.empty:
        print("[WARNING] All candidates were rejected due to stability issues.")
        sys.exit(1)

    # 4. BLAST Validation
    final_candidates = []
    if not args.skip_blast:
        print("\n>>> Step 4: Running NCBI BLAST Validation...")
        final_candidates = run_blast_validation(sorted_df)
    else:
        print("\n>>> Step 4: BLAST Validation Skipped.")
        # Convert DataFrame rows to list of dicts for consistency if skipping blast
        final_candidates = sorted_df.to_dict('records')

    # 5. Export
    if final_candidates:
        print(f"\n>>> Step 5: Exporting {len(final_candidates)} candidates...")
        # If final_candidates is a list of dicts (from blast or conversion), convert to df.
        # But export_candidate_results expects list of dicts or df?
        # Looking at pipeline.py, run_blast_validation returns list of series/dicts.
        # Let's ensure consistency.
        
        # Ensure final_candidates is a DataFrame or list of dicts appropriate for export
        if isinstance(final_candidates, list):
             # If elements are Series, convert to dicts
            final_candidates = [item.to_dict() if isinstance(item, pd.Series) else item for item in final_candidates]
            
        final_df = export_candidate_results(final_candidates)
        
        # Override output path if export_candidate_results doesn't take it (it seems it has internal path logic)
        # We should save it to args.output
        if final_df is not None:
            try:
                final_df.to_excel(args.output, index=False)
                print(f"[SUCCESS] Results saved to {args.output}")
            except Exception as e:
                print(f"[ERROR] Failed to save to {args.output}: {e}")
    else:
        print("[WARNING] No candidates survived validation.")

if __name__ == "__main__":
    main()
