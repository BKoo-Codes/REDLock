from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN
 
# ===== Utility Functions =====
def calculate_gc(seq):
    """
    Calculates the GC content percentage of a DNA sequence.
    """
    gc_count = seq.count("G") + seq.count("C")
    return 100 * gc_count / len(seq) if seq else 0


def compute_oligo_conc(strand1, strand2):
    """
    Computes the effective oligonucleotide concentration for Tm calculation.
    """
    return strand1 - (strand2 / 2) if strand1 >= strand2 else (strand1 + strand2) / 4


def tm_nn(seq, strand1=250, strand2=250, Na=20, Tris=50, Mg=10, dNTPs=1.0, saltcorr=7):
    """
    Calculates the melting temperature (Tm) using Nearest-Neighbor thermodynamics.

    This function utilizes Biopython's Tm_NN engine with custom salt and 
    buffer corrections suitable for RCA/Padlock applications.

    Args:
        seq (str): The DNA sequence to calculate Tm for.
        strand1 (float): Concentration of the first strand (nM). Defaults to 250.
        strand2 (float): Concentration of the second strand (nM). Defaults to 250.
        Na (float): Sodium ion concentration (mM). Defaults to 20.
        Tris (float): Tris buffer concentration (mM). Defaults to 50.
        Mg (float): Magnesium ion concentration (mM). Defaults to 10.
        dNTPs (float): dNTPs concentration (mM). Defaults to 1.0.
        saltcorr (int): Salt correction method (Schildkraut=1, Wetmur=2, SantaLucia=3, 
                        etc.). Defaults to 7 (Owczarzy et al. 2008).

    Returns:
        float: The calculated melting temperature in Celsius (°C). 
               Returns 0 if an error occurs during calculation.
    """
    oligo_conc = compute_oligo_conc(strand1, strand2)
    try:
        return round(Tm_NN(seq, Na=Na, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                           dnac1=oligo_conc, dnac2=oligo_conc, saltcorr=saltcorr), 2)
    except Exception:
        return 0
    
def calc_dg_nn(seq):
    """
    Calculates the Gibbs free energy (delta G) using the Nearest-Neighbor (NN) model.
    """
    # NN thermodynamic parameters for delta G (kcal/mol)
    NN_DG = {
        'AA': -1.00, 'TT': -1.00, 'AT': -0.88, 'TA': -0.58,
        'CA': -1.45, 'TG': -1.45, 'GT': -1.44, 'AC': -1.44,
        'CT': -1.28, 'AG': -1.28, 'GA': -1.30, 'TC': -1.30,
        'CG': -2.17, 'GC': -2.24, 'GG': -1.84, 'CC': -1.84
    }
    seq = seq.upper()
    return round(sum(NN_DG.get(seq[i:i+2], -1.0) for i in range(len(seq)-1)), 2)


def rc(seq):
    """Returns the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


# def is_perfect_complement(seq1, seq2):
#     """
#     Checks if two DNA sequences are perfectly complementary to each other.

#     Args:
#         seq1 (str): The first DNA sequence.
#         seq2 (str): The second DNA sequence.

#     Returns:
#         bool: True if sequences are perfect complements, False otherwise.
#     """
#     return (seq1 == rc(seq2) or seq2 == rc(seq1)) and bool(seq1) and bool(seq2)


# def extract_region_bidir(full_seq, site_seq):
#     """
#     Searches for a specific site within a sequence in both orientations.

#     This function attempts to find the target sequence in its direct form 
#     first, then searches for its reverse complement if the direct form 
#     is not found.

#     Args:
#         full_seq (str): The template sequence to search within.
#         site_seq (str): The target sequence to locate.

#     Returns:
#         tuple[str, str]: A tuple containing (matched_sequence, orientation).
#                          Orientation is either 'direct' or 'rc'.
#                          Returns ('', '') if no match is found.
#     """
#     idx = full_seq.find(site_seq)
#     if idx != -1:
#         return (full_seq[idx:idx+len(site_seq)], "direct")
    
#     rc_site = rc(site_seq)
#     idx_rc = full_seq.find(rc_site)
#     if idx_rc != -1:
#         return (full_seq[idx_rc:idx_rc+len(rc_site)], "rc")
    
#     return ('', '')