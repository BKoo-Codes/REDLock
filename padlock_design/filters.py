from .utils import rc, tm_nn
import re


# ===== Filter Functions =====

def has_long_homopolymer(seq, n=3):
    return any(base * n in seq.upper() for base in "ACGT")


def has_g_quadruplex(seq, g_run=3, loop_max=7):
    """
    True if a canonical G-quadruplex motif is present: >=4 runs of >=g_run G's,
    separated by 1..loop_max nt loops. Also flags the C-rich complement, because
    C-tracts on the padlock become G-tracts on the RCA product. G4/C4 motifs stall
    phi29 polymerase and can impair Cas12a activity.
    """
    seq = seq.upper()
    g_pat = r"G{%d,}(?:.{1,%d}G{%d,}){3,}" % (g_run, loop_max, g_run)
    c_pat = r"C{%d,}(?:.{1,%d}C{%d,}){3,}" % (g_run, loop_max, g_run)
    return bool(re.search(g_pat, seq) or re.search(c_pat, seq))

def has_self_complementarity(seq, min_match=4):
    seq = seq.upper()
    return any(seq[i:i+min_match] in rc(seq) for i in range(len(seq) - min_match + 1))

def has_self_complementarity_global(seq, min_match=6):
    seq_rc = rc(seq)
    for i in range(len(seq) - min_match + 1):
        win = seq[i:i+min_match]
        for j in range(len(seq) - min_match + 1):
            if abs(i - j) < min_match:
                continue
            if win == rc(seq[j:j+min_match]):
                return True
        if win in seq_rc:
            return True
    return False

def has_cross_hybridization(seq1, seq2, min_match=4):
    return any(seq1[i:i+min_match] in rc(seq2) for i in range(len(seq1) - min_match + 1))

def has_tandem_repeat(seq, cutoff=6):
    seq = seq.upper()
    for unit in [2, 3, 4, 5]:
        for i in range(len(seq) - unit + 1):
            motif = seq[i:i+unit]
            # count how many times motif repeats consecutively from i
            reps = 1
            while seq[i + unit*reps : i + unit*(reps+1)] == motif:
                reps += 1
            if unit * reps >= cutoff:
                return True
    return False

def has_gc_block(seq, length=4):
    seq = seq.upper()
    return any(all(base in "GC" for base in seq[i:i+length])
               for i in range(len(seq) - length + 1))

def has_palindrome(seq, cutoff=6):
    seq = seq.upper()
    for length in range(cutoff, len(seq) + 1, 2):   # even lengths only
        for i in range(len(seq) - length + 1):
            window = seq[i:i+length]
            if window == rc(window):
                return True
    return False

# ===== Filter Definitions =====

def is_balanced_base_win(seq, min_pct, max_pct):
    L = len(seq)
    freqs = [seq.count(b) / L for b in "ACGT"]
    return all(min_pct <= f <= max_pct for f in freqs)

def spacer_filter(seq):
    return (
        not has_long_homopolymer(seq, 4) and
        not has_gc_block(seq, 4) and
        not has_self_complementarity(seq) and
        is_balanced_base_win(seq, 0.1, 0.50)
    )

def cas12a_filter(seq, tm_range=(45, 63)):
    """Cas12a core screen. `tm_range` was a hard-coded 45–63 °C window; it is now
    a parameter because the ceiling sits below the median Tm at every supported
    core length, so it silently dominates pool size when the length changes."""
    pam_like = ("TTTA", "TTTC", "TTTG", "TAAA", "CAAA", "GAAA")
    tm_min, tm_max = tm_range
    return (
        tm_min <= tm_nn(seq) <= tm_max and
        not has_long_homopolymer(seq, 4) and
        not has_gc_block(seq, 4) and
        not has_self_complementarity(seq) and
        is_balanced_base_win(seq, 0.1, 0.50) and
        all(pam not in seq for pam in pam_like)
    )


# ===== BLAST Utils =====

def check_contiguous_match(query, subject, min_match=15):
    max_run, run = 0, 0
    for q, s in zip(query, subject):
        if q == s:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run >= min_match

def check_3prime_match(query, subject, min_3p=12):
    return query[-min_3p:] == subject[-min_3p:]

def classify_hit(title):
    if re.search(r"\b(NM_|NR_|XM_|XR_)", title):
        return "RNA"
    elif re.search(r"\b(NC_|NG_|AC_|AP_|AL_|CHR)", title):
        return "Genomic"
    else:
        return "Other"