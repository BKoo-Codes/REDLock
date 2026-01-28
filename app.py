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

# Sidebar - Configuration
with st.sidebar:
    st.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=80)
    st.title("Settings")
    
    st.markdown("### 🎯 Target Parameters")
    target_name = st.text_input("Target miRNA Name", value="hsa-miR-21-5p", help="Name of the target miRNA")
    target_seq = st.text_area("Target miRNA Sequence (RNA/DNA)", value="UAGCUUAUCAGACUGAUGUUGA", height=100, help="Sequence of the target miRNA. U will be automatically converted to T.")

    st.markdown("### ⚙️ Design Options")
    num_candidates = st.slider("Number of Candidates", min_value=1, max_value=50, value=10, help="Number of candidates to generate and validate.")
    dg_threshold = st.number_input("Delta G Threshold (kcal/mol)", value=-3.5, step=0.5, help="Stability threshold. Higher values (less negative) are less stable/better.")
    
    st.markdown("### 🛡️ Validation")
    skip_blast = st.checkbox("Skip NCBI BLAST (Faster)", value=False, help="Skip the time-consuming NCBI BLAST step.")
    rnafold_path = st.text_input("RNAfold Path (Optional)", value="", help="Full path to RNAfold executable if not in system PATH.")

    st.markdown("---")
    run_btn = st.button("🚀 Run Design", type="primary")

# Main Content
st.title("🧬 PadDx_miR Padlock Design Tool")
st.markdown("### Automated design pipeline for RCA-based CRISPR-Cas12a detection systems.")

if "results" not in st.session_state:
    st.session_state.results = None

if run_btn:
    if not target_seq:
        st.error("Please enter a target sequence.")
    else:
        with st.status("Running Design Pipeline...", expanded=True) as status:
            # Step 1: Pools
            st.write("Generating sequence pools...")
            spacer_length = [16, 15, 17]
            cas12a_length = [20, 21, 22]
            
            # Use st.cache_data to speed up pool generation if possible, 
            # but for now we run it directly or rely on pipeline's pickle logic.
            # We need to suppress print outputs from the pipeline or redirect them.
            try:
                spacer_pools, cas12a_pools = generate_pool(spacer_length, cas12a_length)
                st.success("✅ Sequence pools initialized.")
            except Exception as e:
                st.error(f"Failed to generate pools: {e}")
                status.update(label="Error in Step 1", state="error")
                st.stop()

            # Step 2: Assembly
            st.write("Assembling candidates...")
            try:
                candidates = assemble_and_filter_padlocks(
                    target_name, 
                    target_seq, 
                    spacer_pools, 
                    cas12a_pools, 
                    num_to_generate=num_candidates
                )
                
                if not candidates:
                    st.error("No valid candidates generated. Try adjusting filters or pool sizes.")
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
                progress_bar = st.progress(0)
                
                # We need to adapt run_blast_validation to work with streamlit progress if possible
                # For now, we just run it. Ideally, we would rewrite the loop here to update UI.
                # Let's adapt the loop manually here for better UI experience.
                
                passed_candidates = []
                target_count = num_candidates # Try to validate up to the requested number
                
                total_to_check = len(sorted_df)
                
                for idx, row in sorted_df.iterrows():
                    candidate_id = row.get('Candidate ID', idx + 1)
                    st.write(f"Checking Candidate #{candidate_id}...")
                    
                    # Basic checks
                    is_safe = True
                    reasons = []
                    
                    # Cas12a Check
                    cas_hits = run_blast_validation(pd.DataFrame([row]), target_count=1) # Re-using function hackily or rewrite?
                    # The imported run_blast_validation does the whole loop and printing. 
                    # It's better to call the lower-level functions if we want granular control, 
                    # but `run_blast_validation` in pipeline.py prints to stdout.
                    # Let's just use the function as is for simplicity, capturing stdout is hard in streamlit easily.
                    # Actually, let's trust the user to wait or see the logs in terminal. 
                    # For a better UI, we should refactor pipeline.py to yield results. 
                    
                    # To keep it simple and robust:
                    final_candidates_list = run_blast_validation(sorted_df) # This prints to console
                    break # We just run it once
                
                # Since run_blast_validation returns a list of approved candidates (dicts/series)
                # But wait, run_blast_validation in pipeline.py returns `passed_candidates` list.
                # And it handles the loop.
                
                # Re-importing to ensure we use the function correctly
                # In pipeline.py: def run_blast_validation(sorted_results, target_count=5): return passed_candidates
                
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
            st.code(f"Padlock: {padlock_seq}")
            st.text(f"Spacer: {row.get('Spacer (Pad)', 'N/A')} | Cas12a Core: {row.get('Cas12a Core (Pad)', 'N/A')}")
            st.markdown("---")
