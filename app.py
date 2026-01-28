import streamlit as st
import pandas as pd
import os
import io
import time
from padlock_design.pipeline import (
    generate_pool, 
    assemble_and_filter_padlocks, 
    analyze_and_sort_padlocks, 
    run_blast_validation,
    export_candidate_results
)
from padlock_design.visualization import plot_structure

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
    .main {
        background-color: #f8f9fa;
    }
    .stButton>button {
        width: 100%;
        background-color: #4CAF50;
        color: white;
        font-weight: bold;
        border-radius: 8px;
        padding: 0.5rem 1rem;
    }
    .stButton>button:hover {
        background-color: #45a049;
        color: white;
    }
    h1 {
        color: #2c3e50;
    }
    h2, h3 {
        color: #34495e;
    }
    .stDataFrame {
        background-color: white;
        border-radius: 10px;
        padding: 10px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
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
    spacer_length_single = st.number_input("Spacer Length (nt)", min_value=0, max_value=30, value=15, help="Length of the spacer sequence to prioritize.")
    cas12a_length_single = st.number_input("Cas12a Core Length (nt)", min_value=20, max_value=24, value=20, help="Length of the Cas12a core sequence to prioritize.")
    
    gc_spacer = st.slider("Spacer GC Content (%)", 0, 100, (35, 60))
    gc_cas = st.slider("Cas12a GC Content (%)", 0, 100, (35, 65))
    
    pool_samples = st.select_slider("Pool Sample Size", options=[10000, 50000, 100000, 200000, 500000, 1000000], value=200000, help="Number of random sequences to generate for initial screening. Higher = better diversity but slower.")

    st.markdown("### ⚙️ General Options")
    num_candidates = st.slider("Number of Candidates", min_value=1, max_value=50, value=10, help="Number of candidates to generate and validate.")
    dg_threshold = st.number_input("Delta G Threshold (kcal/mol)", value=-3.5, step=0.5, help="Stability threshold. Higher values (less negative) are less stable/better.")
    
    # Advanced Settings Expander
    with st.expander("🛠️ Advanced Filters (Expert)", expanded=False):
        st.markdown("**Constraint Scores (Lower = Stricter)**")
        max_homopolymer = st.number_input("Max Homopolymer Length", 3, 10, 4, help="Maximum allowed repetition of a single base (e.g., AAAAA).")
        score_hairpin_global = st.number_input("Global Hairpin Score", 1, 10, 6, help="Self-complementarity score for the full padlock.")
        score_hairpin_module = st.number_input("Module Hairpin Score", 1, 10, 5, help="Self-complementarity score for Spacer/Cas12a parts.")
        
        st.markdown("---")
        st.markdown("**Cross-Hybridization Scores**")
        score_cross_module = st.number_input("Module-Module Cross-Hyb", 1, 10, 5, help="Interaction between Spacer and Cas12a core.")
        score_cross_mirna = st.number_input("miRNA Arm Cross-Hyb", 1, 10, 3, help="Interaction with miRNA binding arms (Critical).")
        score_cross_global = st.number_input("Global Cross-Hyb", 1, 10, 6, help="General non-specific interaction score.")

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

    skip_blast = st.checkbox("Skip NCBI BLAST (Faster)", value=False, help="Skip the time-consuming NCBI BLAST step.")
    
    # ... existing settings ...

    st.markdown("---")
    run_btn = st.button("🚀 Run Design", type="primary")

# Main Content
st.title("🧬 PadDx_miR Padlock Design Tool")
st.markdown("### Automated design pipeline for RCA-based CRISPR-Cas12a detection systems.")

if "results" not in st.session_state:
    st.session_state.results = None

if run_btn:
    if not target_name or not target_seq:
        st.error("Please enter both Target miRNA Name and Sequence.")
    else:
        # Parse inputs
        spacer_length = [int(spacer_length_single)]
        cas12a_length = [int(cas12a_length_single)]
        
        # Construct Validation Params
        custom_validation_params = {
            "has_long_homopolymer(padlock_seq)": max_homopolymer,
            "has_self_complementarity_global(padlock_seq)": score_hairpin_global,
            "has_tandem_repeat(padlock_seq)": 5, # Fixed for now or add option
            
            "has_self_complementarity(spacer)": score_hairpin_module,
            "has_self_complementarity(cas12a_core)": score_hairpin_module,

            "has_cross_hybridization(spacer, cas12a_core)": score_cross_module,
            "has_cross_hybridization(spacer_rc, cas12a_core)": score_cross_module,
            
            "has_cross_hybridization(mirna_left_rc, spacer)": score_cross_mirna,
            "has_cross_hybridization(mirna_right_rc, cas12a_core)": score_cross_mirna,

            "has_cross_hybridization(padlock_seq, spacer)": score_cross_global,
            "has_cross_hybridization(padlock_seq, cas12a_core)": score_cross_global,
            "has_cross_hybridization(padlock_seq, mirna_left_rc)": score_cross_global,
            "has_cross_hybridization(padlock_seq, mirna_right_rc)": score_cross_global,
            
            "has_cross_hybridization(padlock_seq_rc, mirna_seq)": score_cross_global,
        }

        with st.status("Running Design Pipeline...", expanded=True) as status:
            # Step 1: Pools
            st.write(f"Generating sequence pools (Size: {pool_samples}, Spacer GC: {gc_spacer}%, Cas12a GC: {gc_cas}%)...")
            
            try:
                # Force new pool generation if custom params are likely (simple heuristic: always regen for safety in GUI mode)
                # Or use a different filename based on hash. For now, we reuse the function but pass a temp filename or None to avoid overwriting default cache if desired.
                # Actually, overwriting is fine as it reflects user's latest choice.
                spacer_pools, cas12a_pools = generate_pool(
                    spacer_length, 
                    cas12a_length, 
                    pool_file="custom_pools.pkl", # Use a separate file for GUI custom runs
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
            st.write("Assembling candidates with custom filters...")
            try:
                candidates = assemble_and_filter_padlocks(
                    target_name, 
                    target_seq, 
                    spacer_pools, 
                    cas12a_pools, 
                    num_to_generate=num_candidates,
                    validation_params=custom_validation_params
                )
                
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
            if not skip_blast:
                st.write("Running NCBI BLAST validation (this may take a while)...")
                
                # Re-importing to ensure we use the function correctly
                final_candidates_list = run_blast_validation(sorted_df, target_count=num_candidates)
                
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
    with col1:
        st.metric("Candidates Generated", len(df_res))
    with col2:
        avg_tm = df_res['Padlock Tm'].mean() if 'Padlock Tm' in df_res.columns else 0
        st.metric("Avg. Padlock Tm", f"{avg_tm:.2f} °C")
    with col3:
        avg_dg = df_res['Padlock ΔG (kcal/mol)'].mean() if 'Padlock ΔG (kcal/mol)' in df_res.columns else 0
        st.metric("Avg. Delta G", f"{avg_dg:.2f} kcal/mol")

    # DataFrame
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

    # Detailed View Expander
    with st.expander("🔍 Detailed Sequence View"):
        for i, row in df_res.iterrows():
            st.markdown(f"**#{row.get('Candidate ID', i+1)}**")
            padlock_seq = row.get("Padlock Sequence (5'-P)", "N/A")
            structure = row.get("Secondary Structure", "")
            
            col_det1, col_det2 = st.columns([2, 1])
            
            with col_det1:
                st.code(f"Padlock: {padlock_seq}")
                st.text(f"Spacer: {row.get('Spacer (Pad)', 'N/A')} | Cas12a Core: {row.get('Cas12a Core (Pad)', 'N/A')}")
                if structure:
                    st.text(f"Structure: {structure}")
            
            with col_det2:
                if structure and padlock_seq != "N/A":
                    try:
                        # Clean sequence for visualization
                        clean_seq = padlock_seq.replace("[5'-P]-", "")
                        
                        # Call plot_structure without regions
                        img_buf = plot_structure(
                            clean_seq, 
                            structure, 
                            title=f"Candidate #{row.get('Candidate ID', i+1)}"
                        )
                        st.image(img_buf, caption="Predicted 2D Structure", use_container_width=True)
                    except Exception as e:
                        st.error(f"Could not visualize structure: {e}")
            
            st.markdown("---")