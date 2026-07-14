import random
import itertools
import time
from .utils import calculate_gc, rc, tm_nn, calc_dg_nn
from .filters import *

random.seed(10)


# ===== Candidate Generation =====

def generate_biased_dna(length, at_bias=0.25):
    return ''.join(random.choices(
        "ATGC",
        weights=[at_bias, at_bias, (1 - at_bias) / 2, (1 - at_bias) / 2],
        k=length
    ))


def generate_candidates(lengths, num_samples=200000, gc_range=(35, 65), detailed_filter=None):
    raw_pool = []
    min_gc, max_gc = gc_range

    for L in lengths:
        if L <= 10:
            for tup in itertools.product("ATGC", repeat=L):
                s = "".join(tup)
                if min_gc <= calculate_gc(s) <= max_gc:
                    raw_pool.append(s)
        else:
            for _ in range(num_samples):
                s = generate_biased_dna(L)
                if min_gc <= calculate_gc(s) <= max_gc:
                    raw_pool.append(s)

    unique_pool = list(dict.fromkeys(raw_pool))

    if detailed_filter:
        return [seq for seq in unique_pool if detailed_filter(seq)]
    return unique_pool


# ===== miRNA Arm Split =====

def find_best_split(mirna_name, mirna_seq, min_arm_len=9):
    """Splits miRNA into two arms, minimizing GC deviation and Tm imbalance."""
    mirna_seq = mirna_seq.upper().replace("U", "T")
    max_arm_len = len(mirna_seq) - min_arm_len

    best_score = float('inf')
    best_results = {}

    for l_len in range(min_arm_len, max_arm_len + 1):
        left, right = mirna_seq[:l_len], mirna_seq[l_len:]

        gc_l, gc_r = calculate_gc(left), calculate_gc(right)
        tm_l, tm_r = tm_nn(left), tm_nn(right)
        score = abs(gc_l - 50) + abs(gc_r - 50) + abs(tm_l - tm_r)

        if score < best_score:
            best_score = score
            best_results = {
                "mirna_name":     mirna_name,
                "mirna_seq":      mirna_seq,
                "mirna_left":     left,
                "mirna_right":    right,
                "mirna_left_rc":  rc(left),
                "mirna_right_rc": rc(right),
                "tm_l": tm_l, "tm_r": tm_r,
                "gc_l": gc_l, "gc_r": gc_r,
            }

    if not best_results:
        raise ValueError(
            f"No valid split found for {mirna_name} "
            f"(seq len={len(mirna_seq)}, min_arm_len={min_arm_len})"
        )
    return best_results

# ===== Spacer Generation =====

def generate_short_spacer(length, check_seqs):
    for _ in range(2000):
        if length == 3:
            s = ''.join(random.choices("AT", k=3))
        elif length == 4:
            bases = random.choices("AT", k=4)
            random.shuffle(bases)
            s = ''.join(bases)
        else:
            s = ''.join(random.choices("ATGC", weights=[4, 4, 1, 1], k=length))

        if (s.count('A') + s.count('T')) / length < 0.75: continue
        if has_long_homopolymer(s, 4):                     continue
        if has_self_complementarity(s, 5):                 continue
        if "TTT" in s:                                     continue
        if any(seq[i:i+5] in s or seq[i:i+5] in rc(s)
               for seq in check_seqs for i in range(len(seq) - 5 + 1)): continue
        return s
    return None

# ===== Padlock Assembly =====

def assemble_padlock_candidate(spacer_pools, cas12a_pools, split_results,
                               seed, fixed_spacer=None, fixed_cas12a=None):
    # Structure: [5' Arm RC] - [Spacer] - [Cas12a Core] - [Short Spacer] - [3' Arm RC]
    random.seed(seed)

    mirna_left_rc  = rc(split_results['mirna_left'])
    mirna_right_rc = rc(split_results['mirna_right'])

    spacer = fixed_spacer
    if not spacer:
        for length in range(10, 25):
            if spacer_pools.get(length):
                spacer = random.choice(spacer_pools[length])
                break

    cas12a_core = fixed_cas12a
    if not cas12a_core:
        for length in [20, 21, 22]:
            if cas12a_pools.get(length):
                cas12a_core = random.choice(cas12a_pools[length])
                break

    if not spacer or not cas12a_core:
        return None

    avoid_list = [
        split_results['mirna_left'], split_results['mirna_right'],
        mirna_left_rc, mirna_right_rc,
        spacer, rc(spacer),
        cas12a_core, rc(cas12a_core),
        split_results['mirna_seq'], rc(split_results['mirna_seq']),
    ]

    short_spacer = None
    for _ in range(100):
        ss_len = random.choices([4, 3], weights=[1, 1])[0]
        generated = generate_short_spacer(ss_len, avoid_list)
        if generated:
            short_spacer = generated
            break

    if not short_spacer:
        return None

    padlock_seq = mirna_left_rc + spacer + cas12a_core + short_spacer + mirna_right_rc

    return {
        "padlock_seq":    padlock_seq,
        "spacer":         spacer,
        "cas12a_core":    cas12a_core,
        "short_spacer":   short_spacer,
        "mirna_left_rc":  mirna_left_rc,
        "mirna_right_rc": mirna_right_rc,
    }


# ===== Validation =====

def is_safe_candidate(candidate, val_params, mirna_seq):
    """Returns True if candidate passes all thermodynamic/structural rules."""
    from . import filters as _filters

    left_len  = len(candidate['mirna_left_rc'])
    right_len = len(candidate['mirna_right_rc'])

    inner_seq = candidate['padlock_seq'][left_len:-right_len]
    inner_ext = candidate['padlock_seq'][left_len - 3:-(right_len - 3)]

    if "TTT" in inner_seq:
        return False

    context = {
        **candidate,
        "inner_seq":      inner_seq,
        "inner_ext":      inner_ext,
        "padlock_seq_rc": rc(candidate['padlock_seq']),
        "spacer_rc":      rc(candidate['spacer']),
        "mirna_seq":      mirna_seq,
    }

    for rule_key, cutoff in val_params.items():
        func_name, args_raw = rule_key.split("(", 1)
        arg_names = [name.strip() for name in args_raw.rstrip(")").split(",")]

        func_ref = getattr(_filters, func_name, None)
        if not func_ref:
            continue

        args_values = [context.get(name) for name in arg_names]
        if any(v is None for v in args_values):
            continue

        if func_ref(*args_values, cutoff):
            return False

    return True


def validate_padlocks(spacer_pools, cas12a_pools, split_results, val_para,
                      NUM_PADLOCKS=10, seed_offset=0,
                      fixed_spacer=None, fixed_cas12a=None,
                      max_attempts=1_000_000, time_limit_sec=None, progress_cb=None):
    """Assembles and validates padlock candidates until NUM_PADLOCKS is reached."""
    results, seen_pairs = [], set()
    current_seed = seed_offset
    _t0 = time.time()

    for attempts in range(1, max_attempts + 1):
        if time_limit_sec is not None and (time.time() - _t0) > time_limit_sec:
            print(f"[STOP] Time limit {time_limit_sec}s reached at attempt {attempts}.")
            break
        if progress_cb is not None and attempts % 2000 == 0:
            progress_cb(len(results), NUM_PADLOCKS, attempts, time.time() - _t0)
        candidate = assemble_padlock_candidate(
            spacer_pools, cas12a_pools, split_results, current_seed,
            fixed_spacer=fixed_spacer, fixed_cas12a=fixed_cas12a
        )
        current_seed += 1
        if not candidate:
            continue

        pair_key = (candidate['spacer'], candidate['cas12a_core'])
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        if not is_safe_candidate(candidate, val_para, split_results['mirna_seq']):
            continue

        p_seq = candidate['padlock_seq']
        results.append({
            "Target Name":             split_results['mirna_name'],
            "Padlock Sequence (5'-P)": f"[5'-P]-{p_seq}",
            "Padlock Length":          len(p_seq),
            "Padlock Tm":              round(tm_nn(p_seq), 2),
            "Padlock ΔG (kcal/mol)":   calc_dg_nn(p_seq),
            "GC% (Padlock)":           round(calculate_gc(p_seq), 2),
            "Spacer (Pad)":            candidate['spacer'],
            "Spacer Tm":               round(tm_nn(candidate['spacer']), 2),
            "Cas12a Core (Pad)":       candidate['cas12a_core'],
            "Cas12a Tm":               round(tm_nn(candidate['cas12a_core']), 2),
            "Short spacer":            candidate['short_spacer'],
            "miRNA Sequence":          split_results['mirna_seq'],
            "miRNA 5' Arm (Pad)":      candidate['mirna_left_rc'],
            "miRNA 5' Arm len":        len(candidate['mirna_left_rc']),
            "miRNA 3' Arm (Pad)":      candidate['mirna_right_rc'],
            "miRNA 3' Arm len":        len(candidate['mirna_right_rc']),
            "miRNA Left Arm (RC)":     candidate['mirna_left_rc'],
            "miRNA Right Arm (RC)":    candidate['mirna_right_rc'],
        })
        print(f"[Padlock {len(results)}/{NUM_PADLOCKS}] Generated. (Attempts: {attempts})")
        if progress_cb is not None:
            progress_cb(len(results), NUM_PADLOCKS, attempts, time.time() - _t0)

        if len(results) >= NUM_PADLOCKS:
            break

    if len(results) < NUM_PADLOCKS:
        print(f"[WARNING] Target not reached. Collected {len(results)}/{NUM_PADLOCKS}.")

    return results