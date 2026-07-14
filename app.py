import streamlit as st
import pandas as pd
import os
import io
import time
import json
import datetime
from padlock_design.pipeline import (
    generate_pool, 
    assemble_and_filter_padlocks, 
    analyze_and_sort_padlocks, 
    run_blast_validation,
    export_candidate_results
)
from padlock_design.visualization import plot_structure
from padlock_design.reporter import generate_detailed_report

# ... existing imports ...

# Page Config
st.set_page_config(
    page_title="PadDx_miR Designer",
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
    .stApp code{color:#3ef2b3}

    /* ---- Misc accents ---- */
    hr{border-color:var(--line)}
    a{color:#38bdf8}
</style>
""", unsafe_allow_html=True)

# Helper function to parse list inputs
def parse_int_list(input_str, default):
    try:
        return [int(x.strip()) for x in input_str.split(',')]
    except:
        return default

# Sidebar - Configuration
with st.sidebar:
    st.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=80)
    st.title("Settings")
    
    st.markdown("### 🎯 Target Parameters")
    target_name = st.text_input("Target miRNA Name", value="", help="Name of the target miRNA")
    target_seq = st.text_area("Target miRNA Sequence (RNA/DNA)", value="", height=100, help="Sequence of the target miRNA. U will be automatically converted to T.")

    st.markdown("### 🧬 Structure Design")
    spacer_length_single = st.number_input("Spacer Length (nt)", min_value=0, max_value=30, value=10, help="Length of the spacer sequence to prioritize. (Optimized default: 10)")
    cas12a_length_single = st.number_input("Cas12a Core Length (nt)", min_value=20, max_value=24, value=20, help="Length of the Cas12a core sequence to prioritize.")
    
    gc_spacer = st.slider("Spacer GC Content (%)", 0, 100, (35, 60))
    gc_cas = st.slider("Cas12a GC Content (%)", 0, 100, (35, 65))
    
    pool_samples = st.select_slider("Pool Sample Size", options=[10000, 50000, 100000, 200000, 500000, 1000000], value=200000, help="Number of random sequences to generate for initial screening. Higher = better diversity but slower.")

    st.markdown("### ⚙️ General Options")
    num_candidates = st.slider("Number of Candidates", min_value=1, max_value=50, value=10, help="Number of candidates to generate and validate.")
    random_seed = st.number_input("Random Seed", min_value=0, max_value=999999, value=15, help="Fix for reproducible designs. Same seed + same settings = same candidates.")
    dg_threshold = st.number_input("Delta G Threshold (kcal/mol)", value=-3.5, step=0.5, help="Stability threshold. Higher values (less negative) are less stable/better.")
    
    # Advanced Settings Expander
    with st.expander("🛠️ Advanced Filters (Expert)", expanded=False):
        st.caption("Optimized rule set (from validated notebook). Higher = less strict.")
        st.markdown("**Structural (inner region)**")
        hp_inner   = st.number_input("Inner Poly-N (homopolymer)", 2, 10, 3, help="Min single-base run inside the inner region to reject.")
        sc_padlock = st.number_input("Padlock Hairpin (self-comp)", 4, 12, 7, help="Min self-complementary match in the full padlock to reject.")
        tr_inner   = st.number_input("Inner Tandem Repeat", 4, 12, 6, help="Min tandem-repeat length in inner region to reject.")
        pal_inner  = st.number_input("Inner Palindrome", 4, 12, 6, help="Min palindrome length in extended inner region to reject.")
        sc_module  = st.number_input("Module Hairpin (Spacer/Cas)", 3, 10, 5, help="Min self-complementary match within Spacer or Cas12a core.")

        st.markdown("---")
        st.markdown("**Cross-Hybridization**")
        cross_module = st.number_input("Module & Global Cross-Hyb", 3, 10, 5, help="Spacer↔Cas12a and Padlock↔module interaction.")
        cross_arm    = st.number_input("miRNA-Arm Cross-Hyb", 2, 10, 4, help="Arm↔Cas12a and inner↔arm interaction (arm↔spacer fixed at 3).")
        cross_rca    = st.number_input("RCA↔miRNA Cross-Hyb", 3, 12, 7, help="RCA product re-binding the original miRNA.")
        avoid_g4     = st.checkbox("Avoid G-quadruplex (G4)", value=True, help="Reject GGG-N-GGG-N-GGG-N-GGG and its C-rich complement. G4 stalls phi29 & Cas12a.")

import shutil

# ... existing imports ...

# Helper to find RNAfold
def find_rnafold_path():
    # 1. Check system PATH
    path = shutil.which("RNAfold")
    if path:
        return path
        
    # 2. Check common Windows paths
    common_paths = [
        r"C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe",
        r"C:\Program Files\ViennaRNA Package\RNAfold.exe",
        r"C:\ViennaRNA Package\RNAfold.exe",
        r"C:\ViennaRNA\RNAfold.exe"
    ]
    for p in common_paths:
        if os.path.exists(p):
            return p
            
    return ""

# Sidebar - Configuration
with st.sidebar:
    # ... existing settings ...
    
    st.markdown("### 🛡️ Validation & Stability")
    
    # RNAfold Configuration
    detected_path = find_rnafold_path()
    rnafold_path = st.text_input(
        "RNAfold Path", 
        value=detected_path, 
        help="Full path to RNAfold.exe. If empty, dG calculation will be skipped (0.0)."
    )
    
    if st.button("Test RNAfold Path"):
        if rnafold_path and os.path.exists(rnafold_path):
            st.success("✅ RNAfold found!")
        elif not rnafold_path:
            st.warning("⚠️ Path is empty.")
        else:
            st.error("❌ File not found at this path.")

    st.markdown("---")
    st.markdown("### 💥 BLAST Validation")
    blast_mode = st.radio(
        "Select BLAST Mode",
        options=["Online NCBI BLAST", "Local BLAST", "Skip BLAST"],
        index=0, # Default to Online NCBI BLAST
        help="Choose how to validate candidates against off-target sequences."
    )

    blast_exe_path = ""
    blast_db_path = ""

    if blast_mode == "Local BLAST":
        st.info("For local BLAST, you need to have NCBI BLAST+ installed and a local database prepared.")
        blast_exe_path = st.text_input(
            "Local BLAST Executable Path",
            value="",
            help="Full path to the blastn executable (e.g., C:\\blast\\bin\\blastn.exe)"
        )
        blast_db_path = st.text_input(
            "Local BLAST Database Path",
            value="",
            help="Full path to your local BLAST database (e.g., C:\\blast\\db\\nt)"
        )
        
    st.markdown("---")
    run_btn = st.button("🚀 Run Design", type="primary")

# Main Content
st.markdown("""
<div class="pd-hero">
  <h1>🧬 Pad<span class="pd-grad">Dx</span>_miR — Padlock Design Studio</h1>
  <p>Automated design pipeline for RCA-based CRISPR-Cas12a miRNA detection.</p>
</div>
""", unsafe_allow_html=True)

if "results" not in st.session_state:
    st.session_state.results = None

if run_btn:
    if not target_name or not target_seq:
        st.error("Please enter both Target miRNA Name and Sequence.")
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
                _pool_file = f"pools_s{spacer_length[0]}_c{cas12a_length[0]}_g{gc_spacer[0]}-{gc_spacer[1]}_{gc_cas[0]}-{gc_cas[1]}_n{pool_samples}.pkl"
                spacer_pools, cas12a_pools = generate_pool(
                    spacer_length, 
                    cas12a_length, 
                    pool_file=_pool_file,  # param-specific cache: changing settings regenerates
                    num_samples=pool_samples,
                    gc_range_spacer=gc_spacer,
                    gc_range_cas=gc_cas
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
                    _pmsg.markdown(f"후보 **{found}/{target}** · 시도 {attempts:,}회 · {elapsed:.0f}s")
                except Exception:
                    pass
            try:
                candidates = assemble_and_filter_padlocks(
                    target_name, 
                    target_seq, 
                    spacer_pools, 
                    cas12a_pools, 
                    num_to_generate=num_candidates,
                    seed_offset=random_seed,
                    validation_params=custom_validation_params,
                    max_attempts=300_000,
                    time_limit_sec=180,
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
            try:
                # Handle optional path
                r_path = rnafold_path if rnafold_path.strip() else None
                sorted_df = analyze_and_sort_padlocks(df, dg_threshold=dg_threshold, rnafold_path=r_path)
                
                if sorted_df.empty:
                    st.warning("All candidates were rejected due to stability issues.")
                    status.update(label="Rejected at Stability", state="error")
                    st.stop()
                
                st.success(f"✅ Stability analysis complete. {len(sorted_df)} survivors.")
            except Exception as e:
                st.error(f"Stability analysis failed: {e}")
                status.update(label="Error in Step 3", state="error")
                st.stop()

            # Step 4: BLAST
            final_candidates_list = []
            if blast_mode != "Skip BLAST":
                st.write(f"Running {blast_mode} validation (this may take a while)...")
                
                final_candidates_list = run_blast_validation(
                    sorted_df, 
                    target_count=num_candidates, 
                    blast_mode=blast_mode,
                    blast_exe_path=blast_exe_path,
                    blast_db_path=blast_db_path
                )
                
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
                    "tool": "PadDx_miR — Padlock Design Studio",
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
    st.subheader("🎉 Design Results")
    
    # Summary Metrics
    col1, col2, col3 = st.columns(3)
    df_res = st.session_state.results
    # crRNA / Cas12a targeting (RCA product = reverse complement of padlock)
    if "Cas12a Core (Pad)" in df_res.columns:
        _comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "U": "A"}
        _rc = lambda x: "".join(_comp.get(bb, "N") for bb in str(x).upper()[::-1])
        df_res = df_res.copy()
        df_res["crRNA Spacer (RNA)"] = df_res["Cas12a Core (Pad)"].apply(lambda x: str(x).upper().replace("T", "U"))
        df_res["Cas12a Target on RCA (5to3)"] = df_res["Cas12a Core (Pad)"].apply(_rc)
    with col1:
        st.metric("Candidates Generated", len(df_res))
    with col2:
        avg_tm = df_res['Padlock Tm'].mean() if 'Padlock Tm' in df_res.columns else 0
        st.metric("Avg. Padlock Tm", f"{avg_tm:.2f} °C")
    with col3:
        avg_dg = df_res['Padlock ΔG (kcal/mol)'].mean() if 'Padlock ΔG (kcal/mol)' in df_res.columns else 0
        st.metric("Avg. Delta G", f"{avg_dg:.2f} kcal/mol")

    # DataFrame
    st.caption("Ranked by Penalty Score (lower = better). Warning flags an over-structured 5-prime arm or Cas12a zone.")
    st.dataframe(df_res, use_container_width=True)

    # Download Button
    # Excel conversion
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
        df_res.to_excel(writer, index=False, sheet_name='Padlock Candidates')
    
    st.download_button(
        label="📥 Download Results as Excel",
        data=buffer.getvalue(),
        file_name=f"padlock_design_{target_name}.xlsx",
        mime="application/vnd.ms-excel"
    )

    st.caption("crRNA Spacer = Cas12a Core (as RNA); the Cas12a target on the RCA product is its reverse complement. Confirm orientation for your Cas12a variant (e.g., LbCas12a).")

    if st.session_state.get("run_config"):
        st.download_button(
            label="🧾 Download Run Config (JSON)",
            data=json.dumps(st.session_state.run_config, indent=2, ensure_ascii=False),
            file_name=f"run_config_{target_name}.json",
            mime="application/json",
        )

    # Detailed View Expander
    with st.expander("🔍 Detailed Sequence View"):
        for i, row in df_res.iterrows():
            st.markdown(f"**#{row.get('Candidate ID', i+1)}**")
            padlock_seq_full = row.get("Padlock Sequence (5'-P)", "N/A")
            structure = row.get("Secondary Structure", "")
            
            # --- Column for main info and image ---
            col_det1, col_det2 = st.columns([2, 1])
            
            with col_det1:
                st.code(f"Padlock: {padlock_seq_full}")
                st.text(f"Spacer: {row.get('Spacer (Pad)', 'N/A')} | Cas12a Core: {row.get('Cas12a Core (Pad)', 'N/A')}")
                if structure:
                    st.text(f"Structure: {structure}")
            
            with col_det2:
                if structure and padlock_seq_full != "N/A":
                    try:
                        clean_seq = padlock_seq_full.replace("[5'-P]-", "")
                        img_buf = plot_structure(
                            clean_seq, 
                            structure, 
                            title=f"Candidate #{row.get('Candidate ID', i+1)}"
                        )
                        st.image(img_buf, caption="Predicted 2D Structure", use_container_width=True)
                    except Exception as e:
                        st.warning(f"Could not visualize structure: {e}")

            # --- Full Analysis Report ---
            st.markdown("---")
            
            # Generate and display the detailed report
            report_str = generate_detailed_report(row)
            st.code(report_str)
            
            st.markdown("---")