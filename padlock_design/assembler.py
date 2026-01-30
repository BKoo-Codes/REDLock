import random
import itertools
from .utils import calculate_gc, rc, tm_nn, calc_dg_nn
from .filters import *
random.seed(10)  # 재현성 고정
 
# ===== Candidate Generation =====
def generate_biased_dna(length, at_bias=0.25):
    """
    Generates a random DNA sequence with a specified AT/GC bias.

    Args:
        length (int): The desired length of the DNA sequence.
        at_bias (float): The probability weight for A and T bases. 
                         The remaining weight is split between G and C. 
                         Defaults to 0.25 (resulting in 25% A, 25% T, 25% G, 25% C).

    Returns:
        str: A randomly generated DNA sequence.
    """
    return ''.join(random.choices("ATGC", weights=[at_bias, at_bias, (1-at_bias)/2, (1-at_bias)/2], k=length))

def generate_candidates(lengths, num_samples=200000, gc_range=(35,65), detailed_filter=None):
    """
    Generates a pool of DNA sequence candidates based on length and GC content.

    This function uses exhaustive generation for short sequences (<= 10 nt) 
    and random biased sampling for longer sequences.

    Args:
        lengths (list[int]): A list of sequence lengths to generate.
        num_samples (int): Number of random samples to attempt for each length > 10. 
                           Defaults to 200,000.
        gc_range (tuple[float, float]): The allowed (min, max) GC content percentage. 
                                        Defaults to (35, 65).
        detailed_filter (callable, optional): A secondary function that returns True 
                                              if a sequence passes custom rules.

    Returns:
        list[str]: A list of unique DNA sequences that passed all filters.
    """
    raw_pool_list = []
    min_gc, max_gc = gc_range

    for L in lengths:
        if L <= 10:
            for tup in itertools.product("ATGC", repeat=L):
                s = "".join(tup)
                if min_gc <= calculate_gc(s) <= max_gc:
                    raw_pool_list.append(s)
        else:
            for _ in range(num_samples):
                s = generate_biased_dna(L)
                if min_gc <= calculate_gc(s) <= max_gc:
                    raw_pool_list.append(s)

    unique_pool = list(dict.fromkeys(raw_pool_list))

    if detailed_filter:
        return [seq for seq in unique_pool if detailed_filter(seq)]
    
    return unique_pool

# ===== miRNA arm split =====
def find_best_split(mirna_name, mirna_seq, min_arm_len=9):
    """
    Finds the optimal split point for miRNA into two arms by minimizing GC penalty and Tm difference.

    Args:
        mirna_name (str): Name of the target miRNA.
        mirna_seq (str): miRNA sequence (RNA or DNA string).
        min_arm_len (int): Minimum required length for each arm. Defaults to 9.

    Returns:
        dict: Best split results including sequences, RC arms, and thermodynamic data.
    """
    # Standardize sequence: uppercase and replace U with T
    mirna_seq = mirna_seq.upper().replace("U", "T")
    max_arm_len = len(mirna_seq) - min_arm_len
    
    best_score = float('inf')
    best_results = {}

    for l_len in range(min_arm_len, max_arm_len + 1):
        left, right = mirna_seq[:l_len], mirna_seq[l_len:]
        
        # Calculate physical properties
        gc_l, gc_r = calculate_gc(left), calculate_gc(right)
        tm_l, tm_r = tm_nn(left), tm_nn(right)
        
        # Scoring: lower is better. 
        # Penalizes deviation from 50% GC and large Tm differences between arms.
        gc_penalty = abs(gc_l - 50) + abs(gc_r - 50)
        tm_penalty = abs(tm_l - tm_r)
        score = gc_penalty + tm_penalty
        
        if score < best_score:
            best_score = score
            best_results = {
                "mirna_name": mirna_name,
                "mirna_seq" : mirna_seq,
                "mirna_left": left,
                "mirna_right": right,
                "mirna_left_rc": rc(left),
                "mirna_right_rc": rc(right),
                "tm_l": tm_l, "tm_r": tm_r,
                "gc_l": gc_l, "gc_r": gc_r
            }
            
    return best_results

# ===== Generate Spacer =====
def generate_short_spacer(length, check_seqs):
    """
    Generates a short, high-AT content spacer sequence that avoids specific overlaps.

    Args:
        length (int): Desired length of the short spacer (typically 3 or 4 nt).
        check_seqs (list[str]): List of sequences to avoid matching against 
                                 to prevent unintended cross-hybridization.

    Returns:
        str | None: A validated DNA sequence or None if no valid sequence is 
                    found within 2,000 attempts.
    """
    for _ in range(2000):
        # Weighted generation based on length
        if length == 3:
            s = ''.join(random.choices("AT", k=3))
        elif length == 4:
            bases = random.choices("AT", k=4)
            random.shuffle(bases)
            s = ''.join(bases)       
        else:
            s = ''.join(random.choices('ATGC', weights=[4,4,1,1], k=length))

        if (s.count('A')+s.count('T'))/length < 0.75: continue
        if has_long_homopolymer(s, 4): continue
        if has_self_complementarity(s, 5): continue
        if "TTT" in s: continue 
        if any(seq[i:i+5] in s or seq[i:i+5] in rc(s)
               for seq in check_seqs for i in range(len(seq)-5+1)): continue
        return s
    return None

def assemble_padlock_candidate(spacer_pools, cas12a_pools, split_results, seed, fixed_spacer=None, fixed_cas12a=None):
    """
    Randomly selects components from provided pools to assemble a single padlock probe.
    The final structure follows: [5' Arm RC] - [Spacer] - [Cas12a Core] - [Short Spacer] - [3' Arm RC]
    """
    random.seed(seed)
    
    # Reverse complement of split miRNA arms
    mirna_left_rc = rc(split_results['mirna_left'])
    mirna_right_rc = rc(split_results['mirna_right'])

    # [1] Select Spacer and Cas12a Core
    spacer = fixed_spacer
    if not spacer:
        # Priority-based sampling from length pools
        for length in [16, 17, 15, 18]:
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

    # [2] Generate Valid Short Spacer (Up to 100 attempts)
    # Check list to avoid sequence overlap within the padlock
    avoid_list = [
        split_results['mirna_left'], split_results['mirna_right'], 
        mirna_left_rc, mirna_right_rc, spacer, rc(spacer), 
        cas12a_core, rc(cas12a_core), split_results['mirna_seq'], rc(split_results['mirna_seq'])
    ]
    
    short_spacer = None
    for _ in range(100):
        # Weighted choice for length 3 or 4
        ss_len = random.choices([4, 3], weights=[1, 1])[0]
        generated_ss = generate_short_spacer(ss_len, avoid_list)
        if generated_ss:
            short_spacer = generated_ss
            break
            
    if not short_spacer:
        return None

    # [3] Final Concatenation
    # Structure: [5' Arm] - [Spacer] - [Cas12a Core] - [Short Spacer] - [3' Arm]
    padlock_seq = mirna_left_rc + spacer + cas12a_core + short_spacer + mirna_right_rc
    
    return {
        "padlock_seq": padlock_seq,
        "spacer": spacer,
        "cas12a_core": cas12a_core,
        "short_spacer": short_spacer,
        "mirna_left_rc": mirna_left_rc,
        "mirna_right_rc": mirna_right_rc
    }

def is_safe_candidate(candidate, val_params, mirna_seq):
    """
    Validates the assembled sequence against the provided thermodynamic and structural rules.
    """
    # Context dictionary mapping argument names to actual sequences
    context = {
        **candidate,
        "padlock_seq_rc": rc(candidate['padlock_seq']),
        "spacer_rc": rc(candidate['spacer']),
        "mirna_seq": mirna_seq
    }

    for rule_key, cutoff in val_params.items():
        # Parsing "function_name(arg1, arg2)"
        func_name, args_raw = rule_key.split("(", 1)
        arg_names = [name.strip() for name in args_raw[:-1].split(",")]
        
        # Accessing the function reference from the global scope
        func_ref = globals().get(func_name)
        if not func_ref:
            continue
            
        # Retrieving sequences based on context keys
        args_values = [context.get(name) for name in arg_names]
        if any(val is None for val in args_values):
            continue 

        # Execute validation check
        if func_ref(*args_values, cutoff):
            return False # Failed a validation rule
            
    return True # All checks passed

def validate_padlocks(spacer_pools, cas12a_pools, split_results, val_para,
                      NUM_PADLOCKS=10, seed_offset=0, fixed_spacer=None, fixed_cas12a=None):
    """
    Executes the main loop to assemble, validate, and collect padlock probe candidates.

    The process repeats until either the target number of candidates is reached 
    or the maximum number of attempts (1,000,000) is exceeded.

    Args:
        spacer_pools (dict): Dictionary of candidate spacers grouped by length.
        cas12a_pools (dict): Dictionary of Cas12a core sequences grouped by length.
        split_results (dict): Dictionary containing the split miRNA arm data.
        val_para (dict): Validation parameters/rules for the safety check.
        NUM_PADLOCKS (int): Target number of valid candidates to generate. Defaults to 10.
        seed_offset (int): Starting random seed for reproducibility. Defaults to 0.
        fixed_spacer (str, optional): A pre-determined spacer to use.
        fixed_cas12a (str, optional): A pre-determined Cas12a core to use.

    Returns:
        list[dict]: A list of dictionaries, each containing a fully validated 
                    padlock candidate with calculated thermodynamic metrics.
    """
    results, seen_pairs = [], set()
    current_seed = seed_offset
    max_attempts = 1000000
    attempts = 0

    while len(results) < NUM_PADLOCKS and attempts < max_attempts:
        attempts += 1
        
        # 1. Assemble candidate
        candidate = assemble_padlock_candidate(
            spacer_pools, cas12a_pools, split_results, 
            current_seed, fixed_spacer, fixed_cas12a
        )
        current_seed += 1
        if not candidate: continue

        # Avoid redundant pairs
        pair_key = (candidate['spacer'], candidate['cas12a_core'])
        if pair_key in seen_pairs: continue
        seen_pairs.add(pair_key)

        # 2. Validate candidate
        if not is_safe_candidate(candidate, val_para, split_results['mirna_seq']):
            continue

# 3. Calculate metrics and format result
        p_seq = candidate['padlock_seq']
        final_entry = {
            "Target Name": split_results['mirna_name'],
            "Padlock Sequence (5'-P)": f"[5'-P]-{p_seq}",
            "Padlock Length": len(p_seq),
            "Padlock Tm": round(tm_nn(p_seq), 2),
            "Padlock ΔG (kcal/mol)": calc_dg_nn(p_seq),
            "GC% (Padlock)": round(calculate_gc(p_seq), 2),
            
            "Spacer (Pad)": candidate['spacer'],
            "Spacer Tm": round(tm_nn(candidate['spacer']), 2),
            
            "Cas12a Core (Pad)": candidate['cas12a_core'],
            "Cas12a Tm": round(tm_nn(candidate['cas12a_core']), 2),
            
            "Short spacer": candidate['short_spacer'],
            "miRNA Sequence": split_results['mirna_seq'],
            "miRNA 5' Arm (Pad)": candidate['mirna_left_rc'],
            "miRNA 5' Arm len": len(candidate['mirna_left_rc']),
            "miRNA 3' Arm (Pad)": candidate['mirna_right_rc'],
            "miRNA 3' Arm len": len(candidate['mirna_right_rc']),
            
            # Add keys for reporter
            "miRNA Left Arm (RC)": split_results['mirna_left_rc'],
            "miRNA Right Arm (RC)": split_results['mirna_right_rc'],
        }
        
        results.append(final_entry)
        print(f"[Padlock {len(results)} / {NUM_PADLOCKS}] Generated. (Attempts: {attempts})")

    if len(results) < NUM_PADLOCKS:
        print(f"[WARNING] COULD NOT REACH TARGET COUNT. COLLECTED {len(results)} CANDIDATES.")

    return results
