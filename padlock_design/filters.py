from .utils import rc, tm_nn
import re

# ===== Filter Functions =====
def has_long_homopolymer(seq, n=3): return any(base*n in seq.upper() for base in "ACGT")

def has_self_complementarity(seq, min_match=4):
    seq = seq.upper()
    seq_rc = rc(seq)  
    
    return any(seq[i:i+min_match] in seq_rc for i in range(len(seq) - min_match + 1))

def has_self_complementarity_global(seq, min_match=6):
    seq_rc = rc(seq)
    for i in range(len(seq)-min_match+1):
        win = seq[i:i+min_match]
        for j in range(len(seq)-min_match+1):
            if abs(i-j) < min_match: continue
            if win == rc(seq[j:j+min_match]): return True
        if win in seq_rc: return True
    return False

def has_cross_hybridization(seq1, seq2, min_match=4):
    return any(seq1[i:i+min_match] in rc(seq2) for i in range(len(seq1)-min_match+1))

def has_tandem_repeat(seq, n=6):
    return any(seq.count(seq[i:i+n]) > 1 for i in range(len(seq)-n+1))

def has_gc_block(seq, length=4):
    """
    Check if the sequence contains a stretch of `length` consecutive bases
    composed only of G or C (any combination, e.g., GCGG, GCCG, CGGC, etc.)
    """
    seq = seq.upper()
    for i in range(len(seq) - length + 1):
        window = seq[i:i + length]
        if all(base in "GC" for base in window):
            return True
    return False


# ===== Filter Definitions =====
def is_balanced_base_win(seq, min_pct, max_pct):
    L = len(seq)
    freqs = [seq.count(b)/L for b in "ACGT"]
    return all(min_pct <= f <= max_pct for f in freqs)

def spacer_filter(seq):         #수정 필요
    g_ratio = seq.count('G') / len(seq)
    c_ratio = seq.count('C') / len(seq)
    return (
        38 <= tm_nn(seq) <= 53 and
        not has_long_homopolymer(seq, 4) and
        not has_gc_block(seq, 4) and
        not has_self_complementarity(seq) and
        g_ratio >= 0.10 and  
        c_ratio >= 0.10 and 
        is_balanced_base_win(seq, 0.10, 0.45)
    )

def cas12a_filter(seq):
    pam_like = ("TTTA", "TTTC", "TTTG", "AAAT", "AAAC", "AAAG")
    g_ratio = seq.count('G') / len(seq)
    c_ratio = seq.count('C') / len(seq)
    return (
        all(pam not in seq for pam in pam_like) and
        45 <= tm_nn(seq) <= 65 and      #65도로 수정
        not has_long_homopolymer(seq, 4) and
        not has_gc_block(seq, 4) and
        not has_self_complementarity(seq) and
        g_ratio >= 0.10 and  
        c_ratio >= 0.10 and 
        is_balanced_base_win(seq, 0.10, 0.45)
    )

import re

def check_contiguous_match(query, subject, min_match=15):
    """
    Checks for continuous identical nucleotide runs between two sequences.

    Args:
        query (str): The probe or primer sequence.
        subject (str): The target or off-target sequence from BLAST.
        min_match (int): The minimum threshold for consecutive matches. Defaults to 15.

    Returns:
        bool: True if a continuous run of matches meets or exceeds min_match.
    """
    max_run, run = 0, 0
    for q, s in zip(query, subject):
        if q == s:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run >= min_match

def check_3prime_match(query, subject, min_3p=12):
    """
    Verifies if the 3' terminus of the sequences are identical.
    High 3' homology is a major risk factor for non-specific extension.

    Args:
        query (str): The probe sequence.
        subject (str): The potential off-target sequence.
        min_3p (int): Length of the 3' end to compare. Defaults to 12.

    Returns:
        bool: True if the 3' ends match exactly.
    """
    return query[-min_3p:] == subject[-min_3p:]

def classify_hit(title):
    """
    Categorizes BLAST hits based on NCBI Accession prefixes.

    Classification Mapping:
    - RNA: NM_ (mRNA), NR_ (ncRNA), XM_/XR_ (Predicted RNA)
    - Genomic: NC_ (Chromosome), NG_ (Genomic region), AC_/AP_/AL_, CHR
    - Other: Everything else

    Args:
        title (str): The hit description or accession ID.

    Returns:
        str: Category string ('RNA', 'Genomic', or 'Other').
    """
    if re.search(r"\b(NM_|NR_|XM_|XR_)", title):
        return "RNA"
    elif re.search(r"\b(NC_|NG_|AC_|AP_|AL_|CHR)", title):
        return "Genomic"
    else:
        return "Other"