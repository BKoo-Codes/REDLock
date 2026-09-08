import streamlit as st
import pandas as pd
import os
import sys
import random
import io
import time
import json
import datetime
from padlock_design.pipeline import (
    generate_pool, 
    assemble_and_filter_padlocks, 
    analyze_and_sort_padlocks, 
    run_blast_validation,
    export_candidate_results,
    RNAfoldError,
    RNAfoldNotFoundError,
)
from padlock_design.visualization import plot_structure
from padlock_design.reporter import generate_detailed_report

# Page Config
st.set_page_config(
    page_title="REDLock Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Sora:wght@600;700;800&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');

    :root{
        --cyan:#22d3ee; --signal:#3ef2b3; --violet:#9b7bff; --magenta:#ec5fd0; --coral:#ff6f6f;
        --line:rgba(255,255,255,.09);
    }

    /* ---- App background: deep base + drifting colour glows ---- */
    .stApp{
        background:
            radial-gradient(1000px 700px at 12% -8%, rgba(34,211,238,.10), transparent 60%),
            radial-gradient(900px 700px at 92% 2%, rgba(155,123,255,.14), transparent 60%),
            radial-gradient(800px 800px at 65% 112%, rgba(236,95,208,.10), transparent 60%),
            #080b16;
        font-family:'Inter',sans-serif;
        color:#eaf0ff;
    }
    .stApp h1,.stApp h2,.stApp h3,.stApp h4{font-family:'Sora',sans-serif;letter-spacing:-.02em;color:#eaf0ff}

    /* ---- Custom hero banner ---- */
    .pd-hero{
        border:1px solid var(--line);border-radius:18px;padding:22px 26px;margin-bottom:10px;
        background:linear-gradient(180deg,rgba(255,255,255,.05),rgba(255,255,255,.015));
        box-shadow:0 0 55px -22px rgba(62,242,179,.4);
    }
    .pd-hero h1{margin:0;font-size:1.95rem;font-weight:800}
    .pd-hero p{margin:.45rem 0 0;color:#93a0c4;font-size:.98rem}
    .pd-grad{background:linear-gradient(115deg,#5eead4,#38bdf8 40%,#a78bfa 72%,#f0abfc);
        -webkit-background-clip:text;background-clip:text;-webkit-text-fill-color:transparent}

    /* ---- Sidebar ---- */
    [data-testid="stSidebar"]{background:#0b0f20;border-right:1px solid var(--line)}

    /* ---- Buttons (gradient + glow) ---- */
    .stButton>button,.stDownloadButton>button{
        width:100%;border:1px solid rgba(34,211,238,.5);border-radius:10px;font-weight:700;
        color:#06141c;background:linear-gradient(120deg,#2ee6cf,#38bdf8 50%,#9b7bff);
        box-shadow:0 8px 24px -10px rgba(56,189,248,.7);transition:transform .15s,box-shadow .2s;
    }
    .stButton>button:hover,.stDownloadButton>button:hover{
        transform:translateY(-2px);box-shadow:0 12px 32px -8px rgba(56,189,248,.95);color:#06141c;
    }

    /* ---- Inputs ---- */
    .stTextInput input,.stTextArea textarea,.stNumberInput input{
        background:rgba(255,255,255,.04)!important;border:1px solid var(--line)!important;
        color:#eaf0ff!important;border-radius:8px!important;
    }

    /* ---- Metrics ---- */
    [data-testid="stMetric"]{
        background:linear-gradient(180deg,rgba(255,255,255,.055),rgba(255,255,255,.02));
        border:1px solid var(--line);border-radius:14px;padding:14px 16px;
    }
    [data-testid="stMetricValue"]{color:#3ef2b3}

    /* ---- DataFrame ---- */
    .stDataFrame{border:1px solid var(--line);border-radius:12px;overflow:hidden;
        box-shadow:0 10px 30px -18px rgba(0,0,0,.8)}

    /* ---- Expander ---- */
    [data-testid="stExpander"]{border:1px solid var(--line)!important;border-radius:12px!important;
        background:rgba(255,255,255,.02)}

    /* ---- Code blocks ---- */
    .stApp pre,.stApp code{background:rgba(3,6,16,.65)!important}
    .stApp code{color:#e6edff}

    /* ---- Misc accents ---- */
    hr{border-color:var(--line)}
    a{color:#38bdf8}
    /* ---- Sidebar brand + section headers ---- */
    .pd-side-head{display:flex;align-items:center;gap:11px;margin:.1rem 0 .5rem}
    .pd-side-dot{font-size:1.55rem;filter:drop-shadow(0 0 9px rgba(62,242,179,.6))}
    .pd-side-title{font-family:'Sora',sans-serif;font-weight:800;font-size:1.15rem;line-height:1;
        background:linear-gradient(115deg,#5eead4,#38bdf8 55%,#a78bfa);
        -webkit-background-clip:text;background-clip:text;-webkit-text-fill-color:transparent}
    .pd-side-sub{color:#93a0c4;font-size:.72rem;letter-spacing:.08em;text-transform:uppercase;margin-top:3px}
    .pd-sec{font-family:'Sora',sans-serif;font-size:.78rem;font-weight:700;letter-spacing:.09em;
        text-transform:uppercase;color:#9fb0d9;margin:1rem 0 .4rem;padding-top:.6rem;border-top:1px solid var(--line)}
    [data-testid="stSidebar"] .stSlider label,[data-testid="stSidebar"] .stNumberInput label,
    [data-testid="stSidebar"] .stTextInput label,[data-testid="stSidebar"] .stTextArea label{
        font-size:.82rem;color:#c7d2f0;font-weight:500}
    /* ---- Selected-candidate banner + report ---- */
    .pd-pick{display:flex;align-items:center;gap:12px;flex-wrap:wrap;
        background:linear-gradient(180deg,rgba(56,189,248,.14),rgba(155,123,255,.08));
        border:1px solid rgba(56,189,248,.35);border-left:4px solid #38bdf8;border-radius:12px;
        padding:12px 16px;margin:.2rem 0 .8rem;font-size:1.02rem;color:#eaf0ff;
        box-shadow:0 8px 26px -16px rgba(56,189,248,.7)}
    .pd-pick b{color:#7dd3fc}
    .pd-pick-tag{font-family:'JetBrains Mono',monospace;font-size:.8rem;color:#cbd5f5;
        background:rgba(255,255,255,.06);border:1px solid var(--line);border-radius:999px;padding:2px 10px}
    .pd-pick-warn{font-size:.8rem;color:#06141c;background:#ffb454;border-radius:999px;padding:2px 10px;font-weight:700}
    .pd-pick-ok{font-size:.8rem;color:#06141c;background:#3ef2b3;border-radius:999px;padding:2px 10px;font-weight:700}
    .pd-report{background:rgba(3,6,16,.72);border:1px solid var(--line);border-radius:12px;
        padding:14px 18px;overflow-x:auto;font-family:'JetBrains Mono',monospace;font-size:.84rem;
        line-height:1.5;margin:0}
    .pd-report>div{white-space:pre-wrap;word-break:break-word}
    .pd-r-head{color:#38bdf8;font-weight:800;font-size:.95rem;letter-spacing:.02em;
        margin:1.05rem 0 .4rem;padding-bottom:.28rem;border-bottom:1px solid rgba(56,189,248,.28)}
    .pd-report>.pd-r-head:first-child{margin-top:0}
    .pd-r-rule{height:1px;background:var(--line);margin:.55rem 0}
    .pd-r-gap{height:.35rem}
    .pd-r-key{color:#3ef2b3;font-weight:700;margin:.2rem 0 .1rem}
    .pd-r-item{padding-left:12px;border-left:2px solid rgba(56,189,248,.22);margin:1px 0}
    .pd-r-k{color:#8f9dc2}
    .pd-r-v{color:#eaf0ff;margin-left:8px}
    /* ---- Inspect: prompt, selectbox, banner, detail ---- */
    .pd-selprompt{font-family:'Sora',sans-serif;font-weight:700;font-size:1.02rem;color:#7dd3fc;
        margin:.3rem 0 .35rem;display:flex;align-items:center;gap:8px}
    [data-baseweb="select"]>div{border:1.6px solid rgba(56,189,248,.6)!important;
        background:rgba(56,189,248,.09)!important;border-radius:10px!important;
        box-shadow:0 0 0 3px rgba(56,189,248,.10);min-height:50px}
    [data-baseweb="select"]>div>div{font-size:1.14rem!important;font-weight:800!important;color:#eaf0ff!important}
    .pd-pick-lead{color:#93a0c4;font-size:.92rem}
    .pd-pick-id{font-family:'Sora',sans-serif;font-weight:800;font-size:1.1rem;color:#7dd3fc}
    .pd-detail{margin:.1rem 0 .2rem}
    .pd-pad-label{font-size:.72rem;letter-spacing:.08em;text-transform:uppercase;color:#93a0c4;margin-bottom:5px}
    .pd-pad-seq{font-family:'JetBrains Mono',monospace;font-size:1.06rem;line-height:1.45;color:#eaf0ff;
        background:rgba(56,189,248,.08);border:1px solid rgba(56,189,248,.28);border-radius:10px;
        padding:11px 13px;white-space:nowrap;overflow-x:auto;margin-bottom:10px}
    .pd-fields div{margin:3px 0;font-size:.82rem;font-family:'JetBrains Mono',monospace}
    .pd-fields .pd-k{display:inline-block;min-width:104px;color:#8f9dc2}
    .pd-fields .pd-v{color:#c7d2f0;word-break:break-all}
</style>
""", unsafe_allow_html=True)

# Helper function to parse list inputs
def parse_int_list(input_str, default):
    try:
        return [int(x.strip()) for x in input_str.split(',')]
    except:
        return default

# Sidebar - Configuration
import shutil


def find_rnafold_path():
    # 0. Next to the running .exe / script, and inside the PyInstaller bundle,
    #    so a shipped RNAfold.exe is found and penalty scoring works.
    _dirs = []
    if getattr(sys, "frozen", False):
        _dirs.append(os.path.dirname(sys.executable))
        _dirs.append(getattr(sys, "_MEIPASS", ""))
    _dirs.append(os.path.dirname(os.path.abspath(__file__)))
    for _d in _dirs:
        if not _d:
            continue
        for _n in ("RNAfold.exe", "RNAfold"):
            for _c in (os.path.join(_d, _n), os.path.join(_d, "ViennaRNA", _n)):
                if os.path.exists(_c):
                    return _c

    # 1. Check system PATH
    path = shutil.which("RNAfold")
    if path:
        return path
    # 2. Check common Windows install locations
    common_paths = [
        r"C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe",
        r"C:\Program Files\ViennaRNA Package\RNAfold.exe",
        r"C:\ViennaRNA Package\RNAfold.exe",
        r"C:\ViennaRNA\RNAfold.exe",
    ]
    for p in common_paths:
        if os.path.exists(p):
            return p
    return ""


def _sec(label):
    """Render a styled sidebar section header."""
    st.markdown(f"<div class='pd-sec'>{label}</div>", unsafe_allow_html=True)


with st.sidebar:
    st.markdown(
        "<div class='pd-side-head'>"
        "<span class='pd-side-dot'>🧬</span>"
        "<div><div class='pd-side-title'>REDLock</div>"
        "<div class='pd-side-sub'>Design Settings</div></div>"
        "</div>",
        unsafe_allow_html=True,
    )

    # ─────────────────────────── ① Target ───────────────────────────
    _sec("① Target miRNA")
    target_name = st.text_input(
        "Name", value="", placeholder="e.g. hsa-miR-21-5p",
        help="Name of the target miRNA.",
    )
    target_seq = st.text_area(
        "Sequence (RNA/DNA)", value="", height=90,
        placeholder="e.g. UAGCUUAUCAGACUGAUGUUGA",
        help="Target miRNA sequence. U is auto-converted to T.",
    )
    target_seq = "".join(c for c in target_seq.upper() if c in "ACGUTN")
    _clean = "".join(c for c in target_seq.upper() if c in "ACGUT")
    if _clean:
        _gc = round(sum(c in "GC" for c in _clean) / len(_clean) * 100, 1)
        st.caption(f"Length **{len(_clean)} nt**  ·  GC **{_gc}%**")

    # ────────────────────────── ② Structure ─────────────────────────
    _sec("② Structure Design")
    _c1, _c2 = st.columns(2)
    with _c1:
        spacer_length_single = st.number_input(
            "Spacer (nt)", min_value=0, max_value=30, value=10,
            help="Spacer length to prioritize. Optimized default: 10.",
        )
    with _c2:
        cas12a_length_single = st.number_input(
            "Cas12a Core (nt)", min_value=20, max_value=24, value=20,
            help="Cas12a core length to prioritize.",
        )
    gc_spacer = st.slider("Spacer GC (%)", 0, 100, (40, 60))
    gc_cas = st.slider("Cas12a GC (%)", 0, 100, (40, 60))
    pool_samples = st.select_slider(
        "Pool Sample Size",
        options=[10000, 50000, 100000, 200000, 500000, 1000000],
        value=500000,
        help="Random sequences per length for initial screening. Higher = more diversity, slower.",
    )

    # ────────────────────────── ③ Generation ────────────────────────
    _sec("③ Generation")
    num_candidates = st.number_input(
        "Candidates to generate", min_value=1, max_value=2000, value=2000, step=10,
        help="How many validated padlocks to generate before ranking / BLAST. Up to 2000.",
    )
    _g1, _g2 = st.columns(2)
    with _g1:
        gen_time_min = st.number_input(
            "Time limit (min)", min_value=1, max_value=60, value=5,
            help="Max time for the LOCAL generation step. Safe to increase.",
        )
    with _g2:
        random_seed = st.number_input(
            "Random seed", min_value=0, max_value=999999, value=15,
            help="Seeds POOL generation (notebook used 15). Assembly offset is fixed at 0.",
        )
    dg_threshold = st.number_input(
        "ΔG threshold (kcal/mol)", value=-5.0, step=0.5,
        help="Stability cutoff. Higher (less negative) = less folded / better accessibility.",
    )
    with st.expander("Generation budget (advanced)", expanded=False):
        max_attempts = st.number_input(
            "Max attempts", min_value=10000, max_value=50_000_000, value=5_000_000, step=100000,
            help="Upper bound on assembly attempts. Kept high so the time limit is the effective cap.",
        )

    # ─────────────────── Advanced Filters (expander) ─────────────────
    with st.expander("🛠️ Advanced Filters (Expert)", expanded=False):
        st.caption("Optimized rule set (validated). Higher = less strict.")
        st.markdown("**Structural (inner region)**")
        hp_inner = st.number_input("Inner Poly-N (homopolymer)", 2, 10, 3, help="Min single-base run inside the inner region to reject.")
        sc_padlock = st.number_input("Padlock Hairpin (self-comp)", 4, 12, 7, help="Min self-complementary match in the full padlock to reject.")
        tr_inner = st.number_input("Inner Tandem Repeat", 4, 12, 6, help="Min tandem-repeat length in inner region to reject.")
        pal_inner = st.number_input("Inner Palindrome", 4, 12, 6, help="Min palindrome length in extended inner region to reject.")
        sc_module = st.number_input("Module Hairpin (Spacer/Cas)", 3, 10, 5, help="Min self-complementary match within Spacer or Cas12a core.")
        st.markdown("**Cross-Hybridization**")
        cross_module = st.number_input("Module & Global Cross-Hyb", 3, 10, 5, help="Spacer↔Cas12a and Padlock↔module interaction.")
        cross_arm = st.number_input("miRNA-Arm Cross-Hyb", 2, 10, 4, help="Arm↔Cas12a and inner↔arm interaction (arm↔spacer fixed at 3).")
        cross_rca = st.number_input("RCA↔miRNA Cross-Hyb", 3, 12, 7, help="RCA product re-binding the original miRNA.")
        avoid_g4 = st.checkbox("Avoid G-quadruplex (G4)", value=False, help="Reject GGG-N-GGG-N-GGG-N-GGG and its C-rich complement. G4 stalls phi29 & Cas12a.")

    # ────────────────────────── ④ Validation ────────────────────────
    _sec("④ Validation")
    detected_path = find_rnafold_path()
    st.markdown("**Thermodynamics (RNAfold)**")
    if detected_path:
        st.caption("✅ RNAfold auto-detected on this machine.")
    else:
        st.caption("⚠️ RNAfold not found — the run will stop at the stability step. "
                   "Install ViennaRNA or set the path below.")
    rnafold_path = st.text_input(
        "RNAfold Path", value=detected_path,
        help="Full path to RNAfold.exe. Required: folding ΔG drives both the "
             "stability filter and the penalty score, so the run stops if it is missing.",
    )
    if st.button("Test RNAfold Path"):
        if rnafold_path and os.path.exists(rnafold_path):
            st.success("✅ RNAfold found!")
        elif not rnafold_path:
            st.warning("⚠️ Path is empty.")
        else:
            st.error("❌ File not found at this path.")

    blast_mode = st.radio(
        "BLAST Mode",
        options=["Online NCBI BLAST", "Local BLAST", "Skip BLAST"],
        index=2,
        help="How to validate candidates against off-target sequences.",
    )
    blast_target = st.number_input(
        "BLAST Pass Target", min_value=1, max_value=50, value=3,
        help="Stop BLAST once this many candidates pass (top-ranked screened first). e.g. 3 gives 3 verified candidates.",
    )
    blast_exe_path = ""
    blast_db_path = ""
    if blast_mode == "Local BLAST":
        st.info("Local BLAST requires NCBI BLAST+ installed and a local database.")
        blast_exe_path = st.text_input(
            "Local BLAST Executable Path", value="",
            help="Full path to blastn (e.g. C:\\blast\\bin\\blastn.exe).",
        )
        blast_db_path = st.text_input(
            "Local BLAST Database Path", value="",
            help="Full path to your local BLAST database (e.g. C:\\blast\\db\\nt).",
        )

    st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
    run_btn = st.button("🚀 Run Design", type="primary")

# Main Content
st.markdown("""
<div class="pd-hero">
  <h1>🧬 <span class="pd-grad">RED</span>Lock — Padlock Design Studio</h1>
  <p>Automated design pipeline for RCA-based CRISPR-Cas12a miRNA detection.</p>
</div>
""", unsafe_allow_html=True)

if "results" not in st.session_state:
    st.session_state.results = None

if run_btn:
    _rnafold_ok = bool(rnafold_path.strip()) and (
        os.path.exists(rnafold_path.strip()) or bool(shutil.which(rnafold_path.strip()))
    )
    if not target_name or not target_seq:
        st.error("Please enter both Target miRNA Name and Sequence.")
    elif not _rnafold_ok:
        # Pre-flight check: ΔG = 0.0 used to be substituted silently, which let
        # every candidate pass the stability filter with a bogus penalty score.
        st.error(
            "RNAfold not available — cannot compute folding ΔG. Install ViennaRNA "
            "(https://www.tbi.univie.ac.at/RNA/) or enter a valid path in "
            "**④ Validation → RNAfold Path**, then run again."
        )
    else:
        # Parse inputs
        spacer_length = [int(spacer_length_single)]
        cas12a_length = [int(cas12a_length_single)]
        
        # Construct Validation Params (optimized rule set — inner_seq/inner_ext based)
        custom_validation_params = {
            "has_long_homopolymer(inner_seq)": hp_inner,
            "has_self_complementarity(padlock_seq)": sc_padlock,
            "has_tandem_repeat(inner_seq)": tr_inner,
            "has_palindrome(inner_ext)": pal_inner,

            "has_self_complementarity(spacer)": sc_module,
            "has_self_complementarity(cas12a_core)": sc_module,

            "has_cross_hybridization(spacer, cas12a_core)": cross_module,
            "has_cross_hybridization(spacer_rc, cas12a_core)": cross_module,

            "has_cross_hybridization(mirna_left_rc, spacer)": 3,
            "has_cross_hybridization(mirna_right_rc, cas12a_core)": cross_arm,

            "has_cross_hybridization(padlock_seq, spacer)": cross_module,
            "has_cross_hybridization(padlock_seq, cas12a_core)": cross_module,
            "has_cross_hybridization(inner_seq, mirna_left_rc)": cross_arm,
            "has_cross_hybridization(inner_seq, mirna_right_rc)": cross_arm,

            "has_cross_hybridization(padlock_seq_rc, mirna_seq)": cross_rca,
        }
        if avoid_g4:
            custom_validation_params["has_g_quadruplex(padlock_seq)"] = 3

        with st.status("Running Design Pipeline...", expanded=True) as status:
            # Step 1: Pools
            st.write(f"Generating sequence pools (Size: {pool_samples}, Spacer GC: {gc_spacer}%, Cas12a GC: {gc_cas}%)...")
            
            try:
                # Force new pool generation if custom params are likely (simple heuristic: always regen for safety in GUI mode)
                # Or use a different filename based on hash. For now, we reuse the function but pass a temp filename or None to avoid overwriting default cache if desired.
                # Actually, overwriting is fine as it reflects user's latest choice.
                _pool_file = f"pools_s{spacer_length[0]}_c{cas12a_length[0]}_g{gc_spacer[0]}-{gc_spacer[1]}_{gc_cas[0]}-{gc_cas[1]}_n{pool_samples}_seed{random_seed}.pkl"
                # Fix the global RNG BEFORE pool generation — exactly like the
                # notebook's random.seed(15). The UI seed drives the POOL.
                # generate_pool() re-seeds with pool_seed and records it in the
                # cache fingerprint, so a stale pickle can no longer be reused.
                random.seed(random_seed)
                spacer_pools, cas12a_pools = generate_pool(
                    spacer_length, 
                    cas12a_length, 
                    pool_file=_pool_file,  # param-specific cache: changing settings regenerates
                    num_samples=pool_samples,
                    gc_range_spacer=gc_spacer,
                    gc_range_cas=gc_cas,
                    pool_seed=int(random_seed),
                )
                st.success("✅ Sequence pools initialized.")
            except Exception as e:
                st.error(f"Failed to generate pools: {e}")
                status.update(label="Error in Step 1", state="error")
                st.stop()

            # Step 2: Assembly
            st.write("Assembling candidates with the optimized rule set...")
            _pbar = st.progress(0.0)
            _pmsg = st.empty()
            def _progress(found, target, attempts, elapsed):
                try:
                    _pbar.progress(min(found / target, 1.0) if target else 0.0)
                    _pmsg.markdown(f"Candidates **{found}/{target}** · {attempts:,} attempts · {elapsed:.0f}s")
                except Exception:
                    pass
            try:
                candidates = assemble_and_filter_padlocks(
                    target_name, 
                    target_seq, 
                    spacer_pools, 
                    cas12a_pools, 
                    num_to_generate=num_candidates,
                    seed_offset=0,  # notebook used 0; UI seed drives the POOL, not this
                    validation_params=custom_validation_params,
                    max_attempts=int(max_attempts),
                    time_limit_sec=int(gen_time_min * 60),
                    progress_cb=_progress
                )
                _pbar.progress(1.0)
                
                if not candidates:
                    st.error("No valid candidates generated. Try relaxing filters (increase scores) or increasing pool size.")
                    status.update(label="Failed at Assembly", state="error")
                    st.stop()
                
                df = pd.DataFrame(candidates)
                st.success(f"✅ {len(df)} candidates assembled.")
            except Exception as e:
                st.error(f"Assembly failed: {e}")
                status.update(label="Error in Step 2", state="error")
                st.stop()

            # Step 3: Stability
            st.write("Analyzing thermodynamic stability (RNAfold)...")
            _sbar = st.progress(0.0)
            _smsg = st.empty()
            def _stab_progress(done, total):
                try:
                    _sbar.progress(min(done / total, 1.0) if total else 0.0)
                    _smsg.markdown(f"RNAfold analysis **{done}/{total}**")
                except Exception:
                    pass
            try:
                r_path = rnafold_path.strip() or None
                sorted_df = analyze_and_sort_padlocks(df, dg_threshold=dg_threshold, rnafold_path=r_path, progress_cb=_stab_progress)
                _sbar.progress(1.0)
                
                if sorted_df.empty:
                    st.warning("All candidates were rejected due to stability issues.")
                    status.update(label="Rejected at Stability", state="error")
                    st.stop()
                
                st.success(f"✅ Stability analysis complete. {len(sorted_df)} survivors.")
            except RNAfoldError as e:
                st.error(f"Folding ΔG could not be computed: {e}")
                st.info("Penalty scoring depends on a real ΔG, so the run was stopped "
                        "instead of scoring every candidate as ΔG = 0.")
                status.update(label="Error in Step 3 (RNAfold)", state="error")
                st.stop()
            except Exception as e:
                st.error(f"Stability analysis failed: {e}")
                status.update(label="Error in Step 3", state="error")
                st.stop()

            # Step 4: BLAST
            final_candidates_list = []
            if blast_mode != "Skip BLAST":
                st.write(f"Running {blast_mode} validation — screening top-ranked candidates until {blast_target} pass...")
                if blast_mode == "Online NCBI BLAST":
                    st.caption("⏳ Online NCBI BLAST is a remote query and can take ~30–90 s per module "
                               "(two per candidate). The status line below updates at each step so you can "
                               "tell it is still running.")
                _bbar = st.progress(0.0)
                _bmsg = st.empty()
                _bstat = st.empty()
                def _blast_progress(screened, total, passed, target):
                    try:
                        frac = min(passed / target, 1.0) if target else (screened / total if total else 0.0)
                        _bbar.progress(frac)
                        _bmsg.markdown(f"BLAST screening **{screened}/{total}** · passed **{passed}/{target}**")
                    except Exception:
                        pass
                def _blast_status(msg):
                    try:
                        _bstat.info(f"🔎 {msg}")
                    except Exception:
                        pass
                
                final_candidates_list = run_blast_validation(
                    sorted_df, 
                    target_count=blast_target, 
                    blast_mode=blast_mode,
                    blast_exe_path=blast_exe_path,
                    blast_db_path=blast_db_path,
                    progress_cb=_blast_progress,
                    status_cb=_blast_status
                )
                _bbar.progress(1.0)
                
                if not final_candidates_list:
                    st.warning("No candidates passed BLAST validation.")
                    status.update(label="Failed at BLAST", state="error")
                    st.stop()
                    
                st.success(f"✅ BLAST validation complete. {len(final_candidates_list)} approved.")
            else:
                st.info("BLAST validation skipped.")
                final_candidates_list = sorted_df.to_dict('records')

            # Finalize
            if final_candidates_list:
                # Convert list of dicts/Series to DataFrame
                if isinstance(final_candidates_list[0], pd.Series):
                    final_candidates_list = [c.to_dict() for c in final_candidates_list]
                
                final_df = export_candidate_results(final_candidates_list)
                st.session_state.results = final_df
                st.session_state.run_config = {
                    "tool": "REDLock — Padlock Design Studio",
                    "generated": datetime.datetime.now().isoformat(timespec="seconds"),
                    "target_name": target_name,
                    "target_sequence": target_seq,
                    "random_seed": int(random_seed),
                    "spacer_length": int(spacer_length_single),
                    "cas12a_length": int(cas12a_length_single),
                    "gc_range_spacer": list(gc_spacer),
                    "gc_range_cas": list(gc_cas),
                    "pool_samples": int(pool_samples),
                    "num_candidates": int(num_candidates),
                    "gen_time_min": int(gen_time_min),
                    "max_attempts": int(max_attempts),
                    "blast_target": int(blast_target),
                    "dg_threshold": float(dg_threshold),
                    "avoid_g4": bool(avoid_g4),
                    "blast_mode": blast_mode,
                    "validation_params": dict(custom_validation_params),
                }
                status.update(label="Design Complete!", state="complete", expanded=False)
            else:
                st.error("Something went wrong during finalization.")
                st.stop()

# Results Display
if st.session_state.results is not None:
    st.divider()

    df_res = st.session_state.results
    # crRNA / Cas12a targeting (RCA product = reverse complement of padlock)
    if "Cas12a Core (Pad)" in df_res.columns:
        _comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "U": "A"}
        _rc = lambda x: "".join(_comp.get(bb, "N") for bb in str(x).upper()[::-1])
        df_res = df_res.copy()
        df_res["crRNA Spacer (RNA)"] = df_res["Cas12a Core (Pad)"].apply(lambda x: str(x).upper().replace("T", "U"))
        df_res["Cas12a Target on RCA (5to3)"] = df_res["Cas12a Core (Pad)"].apply(_rc)

    # ---- Results header ----
    _tname = target_name if target_name else "target"
    st.markdown(
        "<div class='pd-hero' style='margin-top:6px'>"
        f"<h1 style='font-size:1.4rem'>🎉 Design Results <span class='pd-grad'>{_tname}</span></h1>"
        f"<p>{len(df_res)} candidate(s) · ranked by penalty score (lower = better)</p>"
        "</div>",
        unsafe_allow_html=True,
    )

    # ---- Summary metrics ----
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Candidates", len(df_res))
    with m2:
        best_pen = df_res["Penalty Score"].min() if "Penalty Score" in df_res.columns else None
        st.metric("Best Penalty", f"{best_pen:.4f}" if best_pen is not None else "—")
    with m3:
        avg_tm = df_res["Padlock Tm"].mean() if "Padlock Tm" in df_res.columns else 0
        st.metric("Avg. Padlock Tm", f"{avg_tm:.1f} °C")
    with m4:
        avg_dg = df_res["Padlock ΔG (kcal/mol)"].mean() if "Padlock ΔG (kcal/mol)" in df_res.columns else 0
        st.metric("Avg. ΔG", f"{avg_dg:.2f}")

    # ---- Ranked table (key columns first, all data retained) ----
    _sec("Ranked Candidates")
    st.caption("Lower **Penalty Score** = better. **Warning** flags an over-structured 5′ arm or Cas12a zone.")
    _preferred = [
        "Candidate ID", "Penalty Score", "Warning", "Spacer GC%",
        "Padlock Tm", "Padlock ΔG (kcal/mol)", "Calculated ΔG",
        "Padlock Sequence (5'-P)", "Spacer (Pad)", "Cas12a Core (Pad)",
        "crRNA Spacer (RNA)", "Cas12a Target on RCA (5to3)",
    ]
    _cols = [c for c in _preferred if c in df_res.columns] + [c for c in df_res.columns if c not in _preferred]
    st.dataframe(
        df_res[_cols],
        use_container_width=True,
        hide_index=True,
        column_config={
            "Penalty Score": st.column_config.NumberColumn("Penalty Score", format="%.4f"),
        },
    )

    # ---- Downloads ----
    _sec("Export")
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        df_res.to_excel(writer, index=False, sheet_name="Padlock Candidates")
    _dc, _dpad = st.columns([1, 2])
    with _dc:
        st.download_button(
            label="📥 Results (Excel)",
            data=buffer.getvalue(),
            file_name=f"padlock_design_{target_name}.xlsx",
            mime="application/vnd.ms-excel",
            use_container_width=True,
        )
        if st.session_state.get("run_config"):
            st.download_button(
                label="🧾 Run Config (JSON)",
                data=json.dumps(st.session_state.run_config, indent=2, ensure_ascii=False),
                file_name=f"run_config_{target_name}.json",
                mime="application/json",
                use_container_width=True,
            )
    st.caption("crRNA Spacer = Cas12a Core (as RNA); the Cas12a target on the RCA product is its reverse complement. Confirm orientation for your Cas12a variant (e.g., LbCas12a).")

    # ---- Inspect a single candidate ----
    _sec("Inspect Candidate")
    st.markdown("<div class='pd-selprompt'>🔎 Select a candidate below to view its full report</div>", unsafe_allow_html=True)
    _ids = [df_res.iloc[k].get("Candidate ID", k + 1) for k in range(len(df_res))]
    def _fmt_pen(v):
        try:
            return f"{float(v):.4f}"
        except (TypeError, ValueError):
            return str(v)
    _sel = st.selectbox(
        "Select a candidate to view details",
        options=list(range(len(df_res))),
        format_func=lambda k: f"#{_ids[k]}  ·  penalty {_fmt_pen(df_res.iloc[k].get('Penalty Score', '—'))}",
        label_visibility="collapsed",
    )
    row = df_res.iloc[_sel]
    padlock_seq_full = row.get("Padlock Sequence (5'-P)", "N/A")
    structure = row.get("Secondary Structure", "")

    import html as _html
    _pen = row.get("Penalty Score", "—")
    _warn = str(row.get("Warning", "") or "")
    _badge = (f"<span class='pd-pick-warn'>⚠ {_warn}</span>"
              if _warn and _warn != "OK" else "<span class='pd-pick-ok'>✓ OK</span>")
    st.markdown(
        "<div class='pd-pick'>"
        "<span class='pd-pick-lead'>Viewing</span>"
        f"<span class='pd-pick-id'>Candidate #{_ids[_sel]}</span>"
        f"<span class='pd-pick-tag'>penalty {_fmt_pen(_pen)}</span>"
        f"{_badge}</div>",
        unsafe_allow_html=True,
    )

    # Padlock sequence — full width, single line (no wrap avoids stray copy spaces)
    _pad_disp = str(padlock_seq_full).replace(" ", "")
    st.markdown(
        "<div class='pd-detail'>"
        "<div class='pd-pad-label'>Padlock (5′-P)</div>"
        f"<div class='pd-pad-seq'>{_html.escape(_pad_disp)}</div>"
        "</div>",
        unsafe_allow_html=True,
    )

    col_det1, col_det2 = st.columns([2, 1])
    with col_det1:
        _spacer = str(row.get("Spacer (Pad)", "N/A"))
        _cas = str(row.get("Cas12a Core (Pad)", "N/A"))
        _crrna = str(row.get("crRNA Spacer (RNA)", "")) if "crRNA Spacer (RNA)" in df_res.columns else ""
        _fields = (
            f"<div><span class='pd-k'>Spacer</span><span class='pd-v'>{_html.escape(_spacer)}</span></div>"
            f"<div><span class='pd-k'>Cas12a Core</span><span class='pd-v'>{_html.escape(_cas)}</span></div>"
        )
        if _crrna:
            _fields += f"<div><span class='pd-k'>crRNA Spacer</span><span class='pd-v'>{_html.escape(_crrna)}</span></div>"
        if structure:
            _fields += f"<div><span class='pd-k'>Structure</span><span class='pd-v'>{_html.escape(str(structure))}</span></div>"
        st.markdown(f"<div class='pd-fields'>{_fields}</div>", unsafe_allow_html=True)

    with col_det2:
        if structure and padlock_seq_full != "N/A":
            try:
                clean_seq = padlock_seq_full.replace("[5'-P]-", "")
                img_buf = plot_structure(
                    clean_seq,
                    structure,
                    title=f"Candidate #{row.get('Candidate ID', _sel + 1)}",
                )
                st.image(img_buf, caption="Predicted 2D Structure", use_container_width=True)
            except Exception as e:
                st.warning(f"Could not visualize structure: {e}")

    # Full analysis report
    st.markdown("---")
    report_str = generate_detailed_report(row)
    import html as _html
    import re as _re
    _rhtml = []
    for _ln in report_str.split("\n"):
        _sline = _ln.strip()
        if not _sline:
            _rhtml.append("<div class='pd-r-gap'></div>")
        elif "=====" in _ln:
            _m = _re.search(r"\[(.*?)\]", _ln)
            _title = _m.group(1).strip() if _m else _sline.replace("=", "").strip()
            _rhtml.append(f"<div class='pd-r-head'>{_html.escape(_title)}</div>")
        elif set(_sline) == {"-"}:
            _rhtml.append("<div class='pd-r-rule'></div>")
        elif _sline.startswith(">"):
            _rhtml.append(f"<div class='pd-r-key'>{_html.escape(_sline.lstrip('> ').strip())}</div>")
        else:
            _txt = _sline[1:].strip() if _sline.startswith("-") else _sline
            _icls = " pd-r-item" if _sline.startswith("-") else ""
            if ":" in _txt:
                _k, _v = _txt.split(":", 1)
                _rhtml.append(
                    f"<div class='pd-r-line{_icls}'>"
                    f"<span class='pd-r-k'>{_html.escape(_k.strip())}</span>"
                    f"<span class='pd-r-v'>{_html.escape(_v.strip())}</span></div>"
                )
            else:
                _rhtml.append(f"<div class='pd-r-line{_icls}'>{_html.escape(_txt)}</div>")
    st.markdown("<div class='pd-report'>" + "".join(_rhtml) + "</div>", unsafe_allow_html=True)
