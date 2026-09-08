import pandas as pd
import pickle
import os
import random
import time
import math
import subprocess # Added for local BLAST
import tempfile   # Added for local BLAST
import io # Added for local BLAST XML parsing
from Bio.Blast import NCBIWWW, NCBIXML
from .filters import (
    spacer_filter,
    cas12a_filter,
    classify_hit,
    check_contiguous_match,
    check_3prime_match,
)
from .assembler import find_best_split, validate_padlocks, generate_candidates
from .reporter import print_detailed_analysis, export_candidate_results, print_split_report
from .thermo import run_rnafold_dna, RNAfoldError, RNAfoldNotFoundError

POOL_CACHE_VERSION = 2


def _pool_fingerprint(spacer_priorities, cas12a_priorities, num_samples,
                      gc_range_spacer, gc_range_cas, pool_seed, cas12a_tm_range):
    """Identifies the exact parameter set a cached pool was built from."""
    return {
        "version":         POOL_CACHE_VERSION,
        "spacer_lengths":  sorted(int(L) for L in spacer_priorities),
        "cas12a_lengths":  sorted(int(L) for L in cas12a_priorities),
        "num_samples":     int(num_samples),
        "gc_range_spacer": [float(v) for v in gc_range_spacer],
        "gc_range_cas":    [float(v) for v in gc_range_cas],
        "pool_seed":       None if pool_seed is None else int(pool_seed),
        "cas12a_tm_range": [float(v) for v in cas12a_tm_range],
    }


def generate_pool(spacer_priorities, cas12a_priorities, pool_file="pools_by_length.pkl",
                  num_samples=500000, gc_range_spacer=(35, 60), gc_range_cas=(35, 65),
                  pool_seed=None, cas12a_tm_range=(45, 63)):
    """
    Generate or load spacer and Cas12a pools.

    The cache file now stores the parameter fingerprint alongside the pools and
    is only reused when every parameter matches. Previously any existing pickle
    was loaded blindly, so changing the length or GC range and re-running (the
    CLI always used the same default filename) silently reused the old pools.

    Args:
        pool_seed (int, optional): If given, seeds the global RNG before pool
            generation and takes part in the cache fingerprint.
        cas12a_tm_range (tuple): Allowed Cas12a core Tm window in °C.
    """
    fingerprint = _pool_fingerprint(spacer_priorities, cas12a_priorities,
                                    num_samples, gc_range_spacer, gc_range_cas,
                                    pool_seed, cas12a_tm_range)

    if os.path.exists(pool_file):
        print(f"[INFO] Found pool cache {pool_file}; verifying parameters...")
        try:
            with open(pool_file, "rb") as f:
                cached = pickle.load(f)
            if (isinstance(cached, tuple) and len(cached) == 3
                    and cached[0] == fingerprint):
                print("[INFO] Cache matches current parameters. Loading pools.")
                return cached[1], cached[2]
            print("[INFO] Cache was built with different parameters. Regenerating...")
        except Exception as e:
            print(f"[WARNING] Failed to load pools: {e}. Regenerating...")

    if pool_seed is not None:
        random.seed(int(pool_seed))

    print("[INFO] Generating new pools...")
    spacer_pools = {L: generate_candidates([L], num_samples=num_samples, gc_range=gc_range_spacer, detailed_filter=spacer_filter) for L in spacer_priorities}
    _cas_filter = lambda seq: cas12a_filter(seq, tuple(cas12a_tm_range))
    cas12a_pools = {L: generate_candidates([L], num_samples=num_samples, gc_range=gc_range_cas, detailed_filter=_cas_filter) for L in cas12a_priorities}

    print("Spacer pools:", {k: len(v) for k, v in spacer_pools.items()}, "→ Total:", sum(len(v) for v in spacer_pools.values()))
    print("Cas12a pools:", {k: len(v) for k, v in cas12a_pools.items()}, "→ Total:", sum(len(v) for v in cas12a_pools.values()))

    try:
        with open(pool_file, "wb") as f:
            pickle.dump((fingerprint, spacer_pools, cas12a_pools), f)
    except Exception as e:
        print(f"[ERROR] Failed to save pools: {e}")

    return spacer_pools, cas12a_pools

def assemble_and_filter_padlocks(name, seq, spacer_pools, cas12a_pools, num_to_generate=30,
                                 seed_offset=0,
                                 validation_params=None, max_attempts=1_000_000,
                                 time_limit_sec=None, progress_cb=None, stats=None):
    """
    Main pipeline: find best split and generate validated padlock candidates.
    Uses the optimized (junior) rule set by default. If validation_params is
    provided, it REPLACES the defaults entirely (no silent merge).
    """
    target_name = name
    target_seq = seq.replace("U", "T")  # Ensure DNA sequence

    split_res = find_best_split(target_name, target_seq)
    print_split_report(split_res)

    # ===== Optimized validation rule set (ported from the validated notebook) =====
    DEFAULT_VALIDATION_PARAMS = {
        # 01: Poly-N repeats (inner region only)
        "has_long_homopolymer(inner_seq)": 3,
        "has_g_quadruplex(padlock_seq)": 3,
        # 02: Padlock self-complementarity (hairpin risk)
        "has_self_complementarity(padlock_seq)": 7,
        # 03: Tandem repeats & palindrome (inner)
        "has_tandem_repeat(inner_seq)": 6,
        "has_palindrome(inner_ext)": 6,
        # 04: Internal structure within modules
        "has_self_complementarity(spacer)": 5,
        "has_self_complementarity(cas12a_core)": 5,
        # 05: Cross-hybridization between modules
        "has_cross_hybridization(spacer, cas12a_core)": 5,
        "has_cross_hybridization(spacer_rc, cas12a_core)": 5,
        # 06: Cross-hybridization with miRNA arms
        # All four arm rules share one cutoff. The left-arm/spacer rule used to
        # be pinned at 3 while its three siblings used 4, so the "arm" setting
        # only governed part of the arm chemistry.
        "has_cross_hybridization(mirna_left_rc, spacer)": 4,
        "has_cross_hybridization(mirna_right_rc, cas12a_core)": 4,
        # 07: Global vs. module cross-reactivity (inner-based)
        "has_cross_hybridization(padlock_seq, spacer)": 5,
        "has_cross_hybridization(padlock_seq, cas12a_core)": 5,
        "has_cross_hybridization(inner_seq, mirna_left_rc)": 4,
        "has_cross_hybridization(inner_seq, mirna_right_rc)": 4,
        # 08: RCA product vs. original target
        "has_cross_hybridization(padlock_seq_rc, mirna_seq)": 7,
    }

    # REPLACE (not merge): custom params fully override defaults
    active_params = validation_params if validation_params else DEFAULT_VALIDATION_PARAMS

    results = validate_padlocks(
        spacer_pools,
        cas12a_pools,
        split_res,
        active_params,
        NUM_PADLOCKS=num_to_generate,
        seed_offset=seed_offset,
        max_attempts=max_attempts,
        time_limit_sec=time_limit_sec,
        progress_cb=progress_cb,
        stats=stats,
    )

    return results

def _penalty_features(row):
    """Accessibility-aware penalty (ported from validated notebook). Lower = better."""
    spacer = str(row.get("Spacer (Pad)", ""))
    cas12a = str(row.get("Cas12a Core (Pad)", ""))
    try:
        arm_len = int(row.get("miRNA 5' Arm len", 0))
    except Exception:
        arm_len = 0
    cas_start = arm_len + len(spacer)
    cas_end = cas_start + len(cas12a)
    gc_pct = (sum(1 for b in spacer if b in "GC") / len(spacer) * 100.0) if spacer else 0.0
    try:
        dg = float(row.get("Calculated ΔG", 0.0))
    except Exception:
        dg = 0.0
    st = str(row.get("Secondary Structure", "")).strip()

    def max_paired_run(segment):
        run = mx = 0; prev = ""
        for c in segment:
            if c in "()":
                run = run + 1 if c == prev else 1
                mx = max(mx, run)
            else:
                run = 0
            prev = c
        return mx

    f_gc_low  = 1.0 if gc_pct < 40 else 0.0
    f_gc_high = 1.0 if gc_pct > 60 else 0.0
    f_dg = 1.0 / (1.0 + math.exp(2.0 * (dg - (-4.0))))
    f_arm = 0.0 if dg <= -5.0 else (1.0 if max_paired_run(st[0:arm_len]) >= 4 else 0.0)
    f_cas = 0.0 if dg <= -5.0 else (1.0 if max_paired_run(st[cas_start:cas_end]) >= 6 else 0.0)

    W = {"gc_low": 4.176, "gc_high": 2.453, "dg": 3.026, "arm": 5.172, "cas": 2.854}
    penalty = (W["gc_low"] * f_gc_low + W["gc_high"] * f_gc_high +
               W["dg"] * f_dg + W["arm"] * f_arm + W["cas"] * f_cas)
    warn = []
    if f_arm: warn.append("arm")
    if f_cas: warn.append("cas")
    return gc_pct, penalty, ("+".join(warn) if warn else "OK")


def calc_penalty_score(row):
    """EXACT copy of the notebook's calc_penalty_score (returns pd.Series, joined below).
    Kept byte-for-byte identical so the sort input DataFrame is built the same way
    as the notebook -> identical quicksort tie order."""
    spacer  = str(row['Spacer (Pad)'])
    cas12a  = str(row['Cas12a Core (Pad)'])

    arm_len   = int(row["miRNA 5' Arm len"])
    cas_start = arm_len + len(spacer)
    cas_end   = cas_start + len(cas12a)

    gc_pct = sum(1 for b in spacer if b in 'GC') / len(spacer) * 100 if spacer else 0

    dg = row['Calculated ΔG']
    st = str(row['Secondary Structure']).strip()

    def _mpr(segment):
        run = mx = 0; prev = ''
        for c in segment:
            if c in '()':
                run = run + 1 if c == prev else 1
                mx = max(mx, run)
            else:
                run = 0
            prev = c
        return mx

    fgl = 1.0 if gc_pct < 40 else 0.0
    fgh = 1.0 if gc_pct > 60 else 0.0
    fdg = 1.0 / (1.0 + math.exp(2.0 * (dg - (-4.0))))
    far = 0.0 if dg <= -5.0 else (1.0 if _mpr(st[0:arm_len]) >= 4 else 0.0)
    fcs = 0.0 if dg <= -5.0 else (1.0 if _mpr(st[cas_start:cas_end]) >= 6 else 0.0)

    W = {'w_gc_low':4.176,'w_gc_high':2.453,'w_dg':3.026,'w_arm':5.172,'w_cas':2.854}
    penalty = (W['w_gc_low']*fgl + W['w_gc_high']*fgh + W['w_dg']*fdg +
               W['w_arm']*far + W['w_cas']*fcs)

    warn = []
    if far == 1: warn.append('arm')
    if fcs == 1: warn.append('cas')

    return pd.Series({
        'Spacer GC%'    : round(gc_pct, 1),
        'f_gc_low'      : fgl,
        'f_gc_high'     : fgh,
        'f_dg'          : fdg,
        'f_arm'         : far,
        'f_cas'         : fcs,
        'Penalty Score' : penalty,
        'Warning'       : '+'.join(warn) if warn else 'OK',
    })


def analyze_and_sort_padlocks(df, dg_threshold=-3.5, rnafold_path=None, progress_cb=None):
    """
    Analyzes sequences and filters out candidates with ΔG <= threshold (too stable).
    Returns a filtered dataframe sorted by stability.
    """
    
    if 'Candidate ID' not in df.columns:
        df.insert(0, 'Candidate ID', range(1, len(df) + 1))
    
    results_struct = []
    results_dg = []

    print(f"\n{'[ DNA Padlock Batch Analysis & Filtering ]':^75}")
    print(f"{'ID':<4} | {'Target Name':<15} | {'ΔG (kcal/mol)':<15} | {'Status'}")
    print("-" * 85)

    total = len(df)
    done = 0
    for idx, row in df.iterrows():
        raw_seq = str(row["Padlock Sequence (5'-P)"])
        clean_seq = raw_seq.replace("[5'-P]-", "").strip()

        # No fallback to ΔG = 0.0 here: 0.0 sits above every sane stability
        # threshold, so a missing or broken RNAfold used to let every candidate
        # PASS with a meaningless penalty score. Stop instead.
        try:
            struct, dg = run_rnafold_dna(clean_seq, rnafold_executable=rnafold_path)
        except RNAfoldNotFoundError:
            raise
        except RNAfoldError as e:
            raise RNAfoldError(
                f"Folding ΔG failed for candidate "
                f"{row.get('Candidate ID', idx)}: {e}"
            ) from e

        dg = float(dg)
        results_struct.append(struct)
        results_dg.append(dg)
        
        status = "PASS" if dg > dg_threshold else "REJECTED (Too Stable)"
        print(f"{row['Candidate ID']:<4} | {row['Target Name']:<15} | {dg:>10.2f}     | {status}")
        done += 1
        if progress_cb is not None:
            progress_cb(done, total)

    df['Calculated ΔG'] = results_dg
    df['Secondary Structure'] = results_struct
    
    initial_count = len(df)
    filtered_df = df[df['Calculated ΔG'] > dg_threshold].copy()
    rejected_count = initial_count - len(filtered_df)
    
    # ===== Penalty scoring & sorting — EXACT notebook reproduction =====
    # Built identically to the notebook (join of a Series-returning apply, then a
    # plain sort_values with no kind=), so the quicksort tie order matches the
    # notebook byte-for-byte instead of diverging (14 vs 27).
    if not filtered_df.empty:
        scored_df = filtered_df.join(filtered_df.apply(calc_penalty_score, axis=1))
        sorted_df = scored_df.sort_values(by='Penalty Score', ascending=True).reset_index(drop=True)
    else:
        # Keep the column schema stable even when nothing survives, so callers can
        # inspect the frame without a KeyError.
        sorted_df = filtered_df.reindex(
            columns=list(filtered_df.columns) + ['Spacer GC%', 'f_gc_low', 'f_gc_high',
                                                 'f_dg', 'f_arm', 'f_cas',
                                                 'Penalty Score', 'Warning'])
    
    print("-" * 85)
    print(f"Filtering complete!")
    print(f" - Total Processed: {initial_count}")
    print(f" - Rejected (ΔG < {dg_threshold}): {rejected_count}")
    print(f" - Survivors: {len(sorted_df)}")
    print("-" * 85 + "\n")
    
    return sorted_df

def _run_local_blastn_and_parse(query_seq, label, blast_exe_path, blast_db_path, expect=100, word_size=7, hitlist_size=20, min_ident=0.90):
    """
    Executes local blastn and parses the XML output.
    """
    if not os.path.exists(blast_exe_path):
        print(f"[ERROR] Local BLAST executable not found at: {blast_exe_path}")
        return None  # unverified -> reject (never pass a candidate we could not screen)
    
    # Check if blast_db_path is a directory or a file. For local blast, it expects
    # the base name of the database, not necessarily the path to a .nhr/.nin/.nsq file.
    # We can do a rudimentary check for existence of one of the files
    db_exists = False
    for ext in ['.nhr', '.nin', '.nsq']: # Common extensions for nucleotide BLAST databases
        if os.path.exists(blast_db_path + ext):
            db_exists = True
            break
    
    if not db_exists:
        print(f"[ERROR] Local BLAST database not found at: {blast_db_path}. Expecting files like {blast_db_path}.nhr, etc.")
        return None  # unverified -> reject (never pass a candidate we could not screen)

    print(f"\n--- [{label}] | DB: {blast_db_path} (Local) ---")
    print("Query sequence:", query_seq)

    hits = {"perfect": {"RNA": [], "Genomic": [], "Other": []},
            "partial": {"RNA": [], "Genomic": [], "Other": []},
            "contig":  {"RNA": [], "Genomic": [], "Other": []},
            "3prime":  {"RNA": [], "Genomic": [], "Other": []}}

    # Use temporary files for query and output
    with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as query_file, \
         tempfile.NamedTemporaryFile(mode='r+', delete=False, encoding='utf-8') as output_file:
        
        query_file.write(f">query\n{query_seq}")
        query_file_name = query_file.name
        output_file_name = output_file.name

    cmd = [
        blast_exe_path,
        "-query", query_file_name,
        "-db", blast_db_path,
        "-outfmt", "5",  # XML output
        "-out", output_file_name,
        "-evalue", str(expect),  # command-line blastn uses -evalue (NOT -expect, which is the Biopython/qblast name)
        "-word_size", str(word_size),
        "-max_target_seqs", str(hitlist_size),
        "-task", "blastn-short" # Good for short sequences like oligos
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        with open(output_file_name, 'r', encoding='utf-8') as f:
            blast_output = f.read()

        if blast_output:
            # Need to handle potential empty BLAST results XML which causes an error with NCBIXML.read
            if "<Iteration>" not in blast_output: # Simple check if there are any hits
                print(f"[INFO] No BLAST hits found for {label} in local database.")
                return hits # Return empty hits if no iterations (no hits)

            blast_record = NCBIXML.read(io.StringIO(blast_output))
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps: # Corrected from blast_record.hsps
                    subj = hsp.sbjct
                    cat = classify_hit(alignment.title)  # same RNA/Genomic/Other classification as online

                    # 1. Perfect Match
                    if hsp.identities == len(query_seq) and hsp.align_length == len(query_seq):
                        hits["perfect"][cat].append(alignment.title)

                    # 2. Partial Match (High Identity >= min_ident, with coverage) -- matches online
                    identity_ratio = hsp.identities / hsp.align_length if hsp.align_length > 0 else 0
                    if min_ident <= identity_ratio < 1.0 and hsp.align_length >= int(min_ident * len(query_seq)):
                        hits["partial"][cat].append(
                            f"{alignment.title} (identity={identity_ratio:.2f}, mismatches={hsp.align_length-hsp.identities})"
                        )

                    # 3. Contiguous Match Check (consecutive >= 15 nt) -- was MISSING in local
                    if check_contiguous_match(hsp.query, subj, min_match=15):
                        hits["contig"][cat].append(alignment.title)

                    # 4. 3' Terminal Match (3' end >= 13 nt) -- was MISSING in local
                    if check_3prime_match(query_seq, subj, min_3p=13):
                        hits["3prime"][cat].append(alignment.title)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Local BLAST command failed: {e.stderr}")
        print("[ERROR] Candidate will be REJECTED (unverified).")
        return None  # unverified -> evaluate_candidate_safety treats as unsafe (reject)
    except Exception as e:
        print(f"[ERROR] Local BLAST execution or parsing failed: {e}")
        print("[ERROR] Candidate will be REJECTED (unverified).")
        return None  # unverified -> reject, never silently pass
    finally:
        # Clean up temporary files
        if os.path.exists(query_file_name):
            os.remove(query_file_name)
        if os.path.exists(output_file_name): # Output file might not be created if subprocess fails early
            os.remove(output_file_name)

    return hits



# NCBI BLAST Functions
# classify_hit / check_contiguous_match / check_3prime_match used to be
# re-implemented here, byte-identical to the versions in filters.py. They are
# now imported from filters.py at the top of this module so the online and
# local BLAST paths can never drift apart.

def run_blastn_and_parse(query_seq, label, db="refseq_rna", expect=100, word_size=7, hitlist_size=20, min_ident=0.90, status_cb=None):
    """NCBI BLAST 실행 후 결과 파싱"""
    print(f"\n--- [{label}] | DB: {db}, txid9606[ORGN] ---")
    print("Query sequence:", query_seq)

    def _say(msg):
        # surface a human-readable status to the UI (best-effort)
        if status_cb is not None:
            try:
                status_cb(msg)
            except Exception:
                pass

    blast_record = None
    max_retries = 3  # 최초 1회 + 재시도 2회
    for attempt in range(1, max_retries + 1):
        t0 = time.time()
        try:
            _say(f"[{label}] contacting NCBI BLAST (attempt {attempt}/{max_retries}) — "
                 f"this usually takes 30–90 s per query, please wait…")
            result_handle = NCBIWWW.qblast("blastn", db, query_seq,
                                           entrez_query="txid9606[ORGN]",
                                           expect=expect,
                                           word_size=word_size,
                                           hitlist_size=hitlist_size)
            _say(f"[{label}] response received after {time.time()-t0:.0f} s — parsing…")
            blast_record = NCBIXML.read(result_handle)
            result_handle.close()
            break  # 성공하면 루프 탈출
        except Exception as e:
            print(f"[WARN] BLAST attempt {attempt}/{max_retries} failed: {e}")
            if attempt < max_retries:
                _say(f"[{label}] network issue after {time.time()-t0:.0f} s "
                     f"(attempt {attempt}/{max_retries}); retrying in 10 s…")
                print("       Try again in 10 seconds...")
                time.sleep(10)
            else:
                _say(f"[{label}] BLAST failed after {max_retries} attempts — candidate rejected (unverified).")
                print(f"[ERROR] BLAST failed after {max_retries} attempts. "
                      f"Candidate will be REJECTED (unverified).")
                # None을 반환 -> evaluate_candidate_safety가 '판정불가 = 탈락'으로 처리
                return None

    hits = {"perfect": {"RNA": [], "Genomic": [], "Other": []},
            "partial": {"RNA": [], "Genomic": [], "Other": []},
            "contig":  {"RNA": [], "Genomic": [], "Other": []},
            "3prime":  {"RNA": [], "Genomic": [], "Other": []}}

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            subj = hsp.sbjct
            cat = classify_hit(alignment.title)

            # 1. Perfect Match Check
            if hsp.identities == len(query_seq) and hsp.align_length == len(query_seq):
                hits["perfect"][cat].append(alignment.title)

            # 2. Partial Match Check (High Identity > 90%)
            identity_ratio = hsp.identities / hsp.align_length if hsp.align_length > 0 else 0
            if min_ident <= identity_ratio < 1.0 and hsp.align_length >= int(min_ident * len(query_seq)):
                hits["partial"][cat].append(
                    f"{alignment.title} (identity={identity_ratio:.2f}, mismatches={hsp.align_length-hsp.identities})"
                )

            # 3. Contiguous Match Check (Consecutive 15nt)
            if check_contiguous_match(hsp.query, subj, min_match=15):
                hits["contig"][cat].append(alignment.title)

            # 4. 3' Terminal Match (3' end 13nt)
            if check_3prime_match(query_seq, subj, min_3p=13):
                hits["3prime"][cat].append(alignment.title)

    return hits

def evaluate_candidate_safety(hits):
    # BLAST가 재시도 끝에 실패하면 hits=None -> 검증 불가 -> 안전하지 않음(탈락)으로 처리
    if hits is None:
        return False, ["BLAST UNVERIFIED (network failure after retries)"]

    reasons = []
    
    # 1. Perfect Match Check
    for cat in hits["perfect"]:
        if hits["perfect"][cat]:
            reasons.append(f"PERFECT MATCH FOUND IN {cat}")

    # 2. Partial Match Check (High Identity > 90%)
    for cat in hits["partial"]:
        if hits["partial"][cat]:
            reasons.append(f"HIGH IDENTITY PARTIAL MATCH (>90%) IN {cat}")

    # 3. Contiguous Match Check (Consecutive 15nt)
    for cat in hits["contig"]:
        if hits["contig"][cat]:
            reasons.append(f"CONTIGUOUS MATCH (>= 15nt) IN {cat}")

    # 4. 3' Terminal Match (3' end 13nt)
    for cat in hits["3prime"]:
        if hits["3prime"][cat]:
            reasons.append(f"Critical: 3' terminal match >= 13nt in {cat}")

    is_safe = len(reasons) == 0
    return is_safe, reasons

def run_blast_validation(sorted_results, target_count=5, blast_mode="Online NCBI BLAST", blast_exe_path="", blast_db_path="", progress_cb=None, status_cb=None):
    passed_candidates = []
    total = len(sorted_results)
    screened = 0

    def _say(msg):
        if status_cb is not None:
            try:
                status_cb(msg)
            except Exception:
                pass

    for idx, row in sorted_results.iterrows():
        candidate_id = row.get('Candidate ID', idx + 1)
        screened += 1
        if progress_cb is not None:
            progress_cb(screened, total, len(passed_candidates), target_count)
        print(f"\n{'='*30} [STARTING VALIDATION: CANDIDATE #{candidate_id}] {'='*30}")

        # Determine which BLAST function to use
        blast_func = None
        if blast_mode == "Online NCBI BLAST":
            blast_func = lambda query_seq, label: run_blastn_and_parse(query_seq, label, status_cb=status_cb)
            print_msg = "Online NCBI BLAST"
        elif blast_mode == "Local BLAST":
            # Validate paths once, up front. The old code re-checked inside the
            # loop and `continue`d on every candidate, so an invalid path looked
            # exactly like "nothing passed BLAST".
            blast_db_base_path = blast_db_path
            if not os.path.exists(blast_exe_path) or \
               (not os.path.exists(blast_db_base_path + '.nhr') and \
                not os.path.exists(blast_db_base_path + '.nin') and \
                not os.path.exists(blast_db_base_path + '.nsq')):
                raise ValueError(
                    f"Local BLAST is not usable. Executable: {blast_exe_path!r}, "
                    f"database base name: {blast_db_base_path!r} (expected "
                    f"{blast_db_base_path}.nhr / .nin / .nsq). Fix the paths or "
                    "choose a different BLAST mode."
                )

            blast_func = lambda query_seq, label: _run_local_blastn_and_parse(
                query_seq, label, blast_exe_path, blast_db_path
            )
            print_msg = "Local BLAST"
        else:
            raise ValueError(
                f"Invalid BLAST mode {blast_mode!r}. Expected 'Online NCBI BLAST', "
                "'Local BLAST', or 'Skip BLAST'."
            )

        # STEP 1: Validate CAS12A Core first
        cas_name = f"CAS12A_CORE_#{candidate_id}"
        cas_seq = row['Cas12a Core (Pad)']

        _say(f"Candidate #{candidate_id} (screened {screened}/{total}, passed {len(passed_candidates)}/{target_count}) "
             f"· STEP 1/2: Cas12a core · {print_msg}…")
        print(f"▶ STEP 1: Running {print_msg} for {cas_name}...")
        _t = time.time()
        cas_hits = blast_func(cas_seq, cas_name)
        _say(f"Candidate #{candidate_id} · STEP 1/2 done in {time.time()-_t:.0f} s")
        is_cas_safe, cas_reasons = evaluate_candidate_safety(cas_hits)

        if not is_cas_safe:
            print(f"\n[REJECTED] CANDIDATE #{candidate_id} FAILED AT CAS12A VALIDATION.")
            print(f"   - REASON(S): {', '.join(cas_reasons)}")
            _say(f"Candidate #{candidate_id} rejected at STEP 1 ({', '.join(cas_reasons)})")
            continue 

        # STEP 2: Validate SPACER Core only if CAS12A passed
        spacer_name = f"SPACER_#{candidate_id}"
        spacer_seq = row['Spacer (Pad)']

        _say(f"Candidate #{candidate_id} · STEP 2/2: Spacer · {print_msg}…")
        print(f"▶ STEP 2: Running {print_msg} for {spacer_name}...")
        _t = time.time()
        spacer_hits = blast_func(spacer_seq, spacer_name)
        _say(f"Candidate #{candidate_id} · STEP 2/2 done in {time.time()-_t:.0f} s")
        is_spacer_safe, spacer_reasons = evaluate_candidate_safety(spacer_hits)

        if not is_spacer_safe:
            print(f"\n[REJECTED] CANDIDATE #{candidate_id} FAILED AT SPACER VALIDATION.")
            print(f"   - REASON(S): {', '.join(spacer_reasons)}")
            continue

        # SUCCESS
        print(f"\n[SUCCESS] CANDIDATE #{candidate_id} IS VERIFIED AS SAFE.")
        passed_candidates.append(row)
        if progress_cb is not None:
            progress_cb(screened, total, len(passed_candidates), target_count)
        
        if len(passed_candidates) >= target_count:
            print(f"\n[INFO] TARGET COUNT ({target_count}) REACHED. TERMINATING VALIDATION PROCESS.")
            break

        if blast_mode == "Online NCBI BLAST":
            time.sleep(2) # Rate limiting for online BLAST

    return passed_candidates
