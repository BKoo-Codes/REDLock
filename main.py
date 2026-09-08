import argparse
import pandas as pd
import os
import sys
from padlock_design.pipeline import (
    generate_pool,
    assemble_and_filter_padlocks,
    analyze_and_sort_padlocks,
    run_blast_validation,
    export_candidate_results,
    RNAfoldError,
)
from padlock_design.thermo import run_rnafold_dna, set_rnafold_path


def main():
    parser = argparse.ArgumentParser(description="REDLock: Padlock Probe Design Tool")

    parser.add_argument("--target_name", type=str, default="hsa-miR-141-3p", help="Target miRNA name")
    parser.add_argument("--target_seq", type=str, default="UAACACUGUCUGGUAAAGAUGG", help="Target miRNA sequence")
    parser.add_argument("--num_candidates", type=int, default=10, help="Number of candidates to assemble before analysis")
    parser.add_argument("--output", type=str, default="results/padlock_results.xlsx", help="Output Excel file path")
    parser.add_argument("--skip_blast", action="store_true", help="Skip NCBI BLAST validation (faster)")
    parser.add_argument("--dg_threshold", type=float, default=-5.0, help="Delta G threshold for stability filtering")
    parser.add_argument("--rnafold_path", type=str, default=None, help="Path to the RNAfold executable")

    parser.add_argument("--spacer_length", type=int, default=10, help="Spacer length to prioritize (notebook default: 10)")
    parser.add_argument("--cas12a_length", type=int, default=20, help="Cas12a core length to prioritize (notebook default: 20)")
    parser.add_argument("--pool_seed", type=int, default=15, help="Seeds POOL generation (notebook used random.seed(15))")
    parser.add_argument("--seed_offset", type=int, default=0, help="Assembly seed start. Notebook used 0. (NOTE: this is NOT the pool seed)")
    parser.add_argument("--pool_samples", type=int, default=500000, help="Random samples per pool length")
    parser.add_argument("--blast_count", type=int, default=5, help="Number of candidates to keep after BLAST validation")
    parser.add_argument("--cas12a_tm_min", type=float, default=45.0, help="Minimum Cas12a core Tm (°C)")
    parser.add_argument("--cas12a_tm_max", type=float, default=63.0, help="Maximum Cas12a core Tm (°C)")
    parser.add_argument("--pool_file", type=str, default=None, help="Pool cache file. Default: a name derived from the pool parameters")

    args = parser.parse_args()

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"\n{'='*50}")
    print(f"REDLock Padlock Design Tool")
    print(f"{'='*50}")
    print(f"Target: {args.target_name}")
    print(f"Sequence: {args.target_seq}")
    print(f"Spacer/Cas12a length: {args.spacer_length} / {args.cas12a_length}")
    print(f"Cas12a Tm window: {args.cas12a_tm_min}-{args.cas12a_tm_max} °C")
    print(f"Pool seed: {args.pool_seed} | Assembly seed offset: {args.seed_offset}")
    print(f"Candidates to assemble: {args.num_candidates}")
    print(f"Output: {args.output}")
    if args.rnafold_path:
        print(f"RNAfold Path: {args.rnafold_path}")
    print(f"{'='*50}\n")

    # 0. Pre-flight: folding ΔG is mandatory (it drives both the stability filter
    #    and the penalty score), so check RNAfold BEFORE spending minutes on
    #    pool generation and assembly.
    print(">>> Step 0: Checking RNAfold availability...")
    try:
        set_rnafold_path(args.rnafold_path)
        _probe_struct, _probe_dg = run_rnafold_dna("ACGTACGTACGTACGT", rnafold_executable=args.rnafold_path)
        print(f"[OK] RNAfold responded (probe ΔG = {_probe_dg:.2f} kcal/mol).")
    except RNAfoldError as e:
        print(f"[ERROR] {e}")
        print("[ERROR] Aborting: candidates cannot be ranked without a real folding ΔG.")
        sys.exit(1)

    # 1. Generate/Load Pools  (single length each, matching the notebook)
    print(">>> Step 1: Initializing Sequence Pools...")
    spacer_length = [args.spacer_length]
    cas12a_length = [args.cas12a_length]
    gc_range_spacer = (40, 60)
    gc_range_cas = (40, 60)
    cas12a_tm_range = (args.cas12a_tm_min, args.cas12a_tm_max)
    # Parameter-specific cache file. The old code always used the default
    # "pools_by_length.pkl", so re-running with a different length silently
    # reused the previous pools. generate_pool() also verifies the fingerprint.
    pool_file = args.pool_file or (
        f"pools_s{args.spacer_length}_c{args.cas12a_length}"
        f"_g{gc_range_spacer[0]}-{gc_range_spacer[1]}"
        f"_{gc_range_cas[0]}-{gc_range_cas[1]}"
        f"_tm{args.cas12a_tm_min:g}-{args.cas12a_tm_max:g}"
        f"_n{args.pool_samples}_seed{args.pool_seed}.pkl"
    )
    spacer_pools, cas12a_pools = generate_pool(
        spacer_length,
        cas12a_length,
        pool_file=pool_file,
        num_samples=args.pool_samples,
        gc_range_spacer=gc_range_spacer,
        gc_range_cas=gc_range_cas,
        pool_seed=args.pool_seed,
        cas12a_tm_range=cas12a_tm_range,
    )

    # 2. Assemble and Filter  (use current signature; default rule set lives in pipeline)
    print("\n>>> Step 2: Assembling and Filtering Candidates...")
    try:
        candidates = assemble_and_filter_padlocks(
            args.target_name,
            args.target_seq,
            spacer_pools,
            cas12a_pools,
            num_to_generate=args.num_candidates,
            seed_offset=args.seed_offset,
        )
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    if not candidates:
        print("[ERROR] No candidates generated. Check parameters.")
        sys.exit(1)

    df = pd.DataFrame(candidates)

    # 3. Analyze, filter by stability, and rank by Penalty Score
    print("\n>>> Step 3: Analyzing Stability & Penalty Scoring...")
    try:
        sorted_df = analyze_and_sort_padlocks(df, dg_threshold=args.dg_threshold, rnafold_path=args.rnafold_path)
    except RNAfoldError as e:
        print(f"[ERROR] {e}")
        print("[ERROR] Stopping: penalty scoring is meaningless without a real folding ΔG.")
        sys.exit(1)

    if sorted_df.empty:
        print("[WARNING] All candidates were rejected due to stability issues.")
        sys.exit(1)

    # 4. BLAST Validation (preserves penalty ranking order from Step 3)
    if not args.skip_blast:
        print("\n>>> Step 4: Running NCBI BLAST Validation...")
        final_candidates = run_blast_validation(sorted_df, target_count=args.blast_count)
    else:
        print("\n>>> Step 4: BLAST Validation Skipped.")
        final_candidates = sorted_df.to_dict("records")

    # 5. Export
    if final_candidates:
        print(f"\n>>> Step 5: Exporting {len(final_candidates)} candidates...")
        if isinstance(final_candidates, list):
            final_candidates = [
                item.to_dict() if isinstance(item, pd.Series) else item
                for item in final_candidates
            ]

        final_df = pd.DataFrame(final_candidates)
        try:
            final_df.to_excel(args.output, index=False)
            print(f"[SUCCESS] Results saved to {args.output}")
        except Exception as e:
            print(f"[ERROR] Failed to save to {args.output}: {e}")
    else:
        print("[WARNING] No candidates survived validation.")


if __name__ == "__main__":
    main()
