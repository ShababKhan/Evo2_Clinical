import streamlit as st
import pandas as pd
import sqlite3
import os
import subprocess # To potentially call external tools
from datetime import datetime
import io # For handling file uploads

# --- Attempt to import optional bioinformatics libraries ---
try:
    import pysam # For reading VCF files
except ImportError:
    st.error("Library `pysam` not found. Please install it: `pip install pysam`")
    st.stop() # Stop execution if pysam is missing

try:
    import requests # For querying web APIs (like GWAS Catalog)
except ImportError:
     st.error("Library `requests` not found. Please install it: `pip install requests`")
     st.stop() # Stop execution if requests is missing


# --- Configuration & Constants ---
DB_DIR = "data/db"
DB_PATH = os.path.join(DB_DIR, "pulmonary_endo_toolkit.db")
GENES_OF_INTEREST = ["GATA2-AS1", "EPAS1", "Other EMT Gene", "PAH Gene", "CTEPH Gene"] # Example list
GWAS_API_URL = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait"


# --- Database Setup (Identical to v1) ---
def init_db():
    """Initializes the SQLite database and creates tables if they don't exist."""
    os.makedirs(DB_DIR, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    # Variants Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        variant_id TEXT UNIQUE NOT NULL, -- rsID or chr:pos:ref:alt
        chromosome TEXT,
        position INTEGER,
        ref_allele TEXT,
        alt_allele TEXT,
        source_dataset TEXT, -- e.g., 'VCF Upload', 'GWAS Catalog'
        related_gene TEXT, -- Can be populated later
        date_added TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """)
    # Scores Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS scores (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        variant_id INTEGER NOT NULL, -- Foreign key to variants.id
        tool TEXT,
        score_type TEXT,
        score_value REAL,
        cell_type_context TEXT,
        date_scored TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (variant_id) REFERENCES variants (id)
    )
    """)
    # Annotations Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS annotations (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        variant_id INTEGER NOT NULL, -- Foreign key to variants.id
        annotation_type TEXT, -- e.g., 'Disease Association', 'Reported Trait', 'Pathway'
        annotation_value TEXT,
        source TEXT, -- e.g., 'GWAS Catalog API', 'Manual Curation'
        details TEXT, -- e.g., p-value, study accession
        date_annotated TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (variant_id) REFERENCES variants (id)
    )
    """)
    conn.commit()
    conn.close()
    print(f"Database initialized at {DB_PATH}")

# --- Database Utilities (Identical to v1) ---
def db_query(query, params=()):
    """Executes a SELECT query and returns results as a DataFrame."""
    try:
        conn = sqlite3.connect(DB_PATH)
        # Use dictionary cursor for easier access by column name if needed later
        # conn.row_factory = sqlite3.Row
        df = pd.read_sql_query(query, conn, params=params)
        conn.close()
        return df
    except sqlite3.Error as e:
        st.error(f"Database query error: {e}")
        return pd.DataFrame()

def db_execute(query, params=()):
    """Executes an INSERT, UPDATE, or DELETE query."""
    conn = None
    try:
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute(query, params)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        # Provide more specific feedback for UNIQUE constraint errors
        if "UNIQUE constraint failed: variants.variant_id" in str(e):
            # Don't show error in UI for duplicates, handle silently or log
            print(f"Skipping duplicate variant: {params[0] if params else 'Unknown'}")
            return None # Indicate skip/failure
        else:
            st.error(f"Database execution error: {e}")
        if conn:
            conn.rollback()
        return None
    finally:
        if conn:
            conn.close()

# --- Functional Data Fetching/Parsing ---

def parse_vcf(uploaded_file):
    """Parses variants from an uploaded VCF file using pysam."""
    st.info(f"Attempting to parse VCF file: {uploaded_file.name}")
    variants_data = []
    try:
        # Read the uploaded file content directly
        # Pysam needs a filename, so we might need to save temporarily or use BytesIO
        # For simplicity here, let's assume pysam can handle the stream-like object
        # Note: Pysam typically works best with file paths. Handling streams might require more complex setup.
        # A more robust way is to save the upload to a temporary file.

        # --- Using a temporary file (more robust) ---
        temp_file_path = os.path.join("/tmp", uploaded_file.name) # Use /tmp or tempfile module
        with open(temp_file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        vcf_file = pysam.VariantFile(temp_file_path)
        # --- End temporary file approach ---

        count = 0
        limit = 1000 # Limit the number of variants parsed in this example
        for record in vcf_file.fetch():
            if count >= limit:
                st.warning(f"Reached parsing limit of {limit} variants from VCF.")
                break
            # Extract relevant information
            # Assuming standard VCF fields. ALT can have multiple alleles.
            # We'll take the first ALT allele for simplicity here.
            alt_alleles = record.alts if record.alts else ['N/A']
            variant_id = record.id if record.id else f"{record.chrom}:{record.pos}:{record.ref}:{alt_alleles[0]}"

            variants_data.append({
                'variant_id': variant_id,
                'chromosome': record.chrom,
                'position': record.pos,
                'ref_allele': record.ref,
                'alt_allele': alt_alleles[0], # Taking the first ALT allele
                'source_dataset': f'VCF Upload: {uploaded_file.name}',
                'related_gene': None # Gene annotation would be a separate step
            })
            count += 1
        vcf_file.close()
        os.remove(temp_file_path) # Clean up temporary file
        st.success(f"Successfully parsed {len(variants_data)} variants from {uploaded_file.name}.")

    except Exception as e:
        st.error(f"Error parsing VCF file '{uploaded_file.name}': {e}")
        # Clean up temp file in case of error
        if 'temp_file_path' in locals() and os.path.exists(temp_file_path):
            os.remove(temp_file_path)
        return pd.DataFrame() # Return empty DataFrame on error

    return pd.DataFrame(variants_data)


def fetch_gwas_catalog_by_trait(trait):
    """Fetches study and variant data from GWAS Catalog API based on reported trait."""
    st.info(f"Querying GWAS Catalog API for trait: {trait}")
    variants_data = []
    annotations_data = []

    params = {'traitName': trait}
    try:
        response = requests.get(GWAS_API_URL, params=params, timeout=30) # Added timeout
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        data = response.json()

        studies = data.get('_embedded', {}).get('studies', [])
        if not studies:
            st.warning(f"No studies found for trait '{trait}' in GWAS Catalog.")
            return pd.DataFrame(), pd.DataFrame()

        st.success(f"Found {len(studies)} studies for trait '{trait}'. Fetching associations...")

        study_count = 0
        variant_limit_per_study = 50 # Limit variants per study for performance

        for study in studies:
            study_count += 1
            if study_count > 10: # Limit number of studies processed
                 st.warning("Reached study processing limit (10).")
                 break

            accession_id = study.get('accessionId')
            associations_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accession_id}/associations"

            try:
                assoc_response = requests.get(associations_url, timeout=30)
                assoc_response.raise_for_status()
                assoc_data = assoc_response.json()
                associations = assoc_data.get('_embedded', {}).get('associations', [])

                variant_count_in_study = 0
                for assoc in associations:
                    if variant_count_in_study >= variant_limit_per_study:
                        break

                    # Extract variant info (can be complex due to structure)
                    snps = assoc.get('loci', [{}])[0].get('strongestRiskAlleles', [])
                    if not snps: continue # Skip if no SNP info

                    snp = snps[0] # Take the first strongest risk allele listed
                    variant_id = snp.get('riskAlleleName') # Often chr:pos-allele
                    rsid = snp.get('rsId')
                    if rsid:
                        variant_id = rsid # Prefer rsID if available

                    # Extract position if possible from riskAlleleName or loci
                    chromosome = None
                    position = None
                    if ':' in variant_id:
                         parts = variant_id.split(':')
                         if len(parts) > 1:
                             chromosome = parts[0]
                             try:
                                 pos_part = parts[1].split('-')[0] # Handle like 'chr1:12345-A'
                                 position = int(pos_part)
                             except (ValueError, IndexError):
                                 position = None # Couldn't parse position

                    # Alleles might need parsing from riskAlleleName or context
                    ref_allele = None # Often not directly available in this part
                    alt_allele = snp.get('riskAlleleName', '').split('-')[-1] if '-' in snp.get('riskAlleleName', '') else None

                    p_value = assoc.get('pValue')

                    variants_data.append({
                        'variant_id': variant_id,
                        'chromosome': chromosome,
                        'position': position,
                        'ref_allele': ref_allele, # May need refinement
                        'alt_allele': alt_allele, # May need refinement
                        'source_dataset': 'GWAS Catalog API',
                        'related_gene': None # Gene info might be elsewhere in API response
                    })

                    # Prepare annotation data
                    annotations_data.append({
                        'variant_id_ref': variant_id, # Use the variant ID for linking later
                        'annotation_type': 'Reported Trait',
                        'annotation_value': trait,
                        'source': f'GWAS Catalog API ({accession_id})',
                        'details': f"p-value: {p_value}" if p_value else None
                    })
                    variant_count_in_study += 1

            except requests.exceptions.RequestException as e_assoc:
                st.warning(f"Could not fetch associations for study {accession_id}: {e_assoc}")
            except Exception as e_proc:
                 st.warning(f"Error processing associations for study {accession_id}: {e_proc}")


    except requests.exceptions.RequestException as e:
        st.error(f"Error querying GWAS Catalog API: {e}")
        return pd.DataFrame(), pd.DataFrame() # Return empty DataFrames on error
    except Exception as e_json:
         st.error(f"Error processing GWAS Catalog response: {e_json}")
         return pd.DataFrame(), pd.DataFrame()

    variants_df = pd.DataFrame(variants_data).drop_duplicates(subset=['variant_id'])
    annotations_df = pd.DataFrame(annotations_data)

    st.success(f"Processed {len(variants_df)} unique variants and {len(annotations_df)} annotations from GWAS Catalog for trait '{trait}'.")
    return variants_df, annotations_df


# --- Placeholder Functions (Unchanged from v1) ---

def run_evo2_scoring(variant_df):
    """Placeholder: Simulates running Evo2 scoring."""
    st.info("Placeholder: Running Evo2 scoring on variants...")
    # In reality: Use subprocess to call Evo2 executable
    scores = pd.DataFrame({
        'variant_id_internal': variant_df['id'], # Use internal DB ID
        'tool': ['Evo2_Placeholder'] * len(variant_df),
        'score_type': ['Functional Impact'] * len(variant_df),
        'score_value': [round(0.1 + i * 0.1, 2) for i in range(len(variant_df))],
        'cell_type_context': ['All'] * len(variant_df)
    })
    st.success("Placeholder: Evo2 scoring complete.")
    return scores

def filter_encode_data(scored_variants_df, cell_type="ENDOS"):
    """Placeholder: Simulates filtering based on ENCODE data."""
    st.info(f"Placeholder: Filtering variants/scores for cell type: {cell_type}")
    # In reality: Query ENCODE API or local data, then filter/update scores
    scored_variants_df['cell_type_context'] = scored_variants_df['cell_type_context'].apply(
        lambda x: cell_type if x == 'All' and hash(str(x)) % 2 == 0 else x
    )
    st.success("Placeholder: ENCODE filtering applied.")
    return scored_variants_df[scored_variants_df['cell_type_context'] == cell_type]

def run_aido_simulation(filtered_variants_df):
    """Placeholder: Simulates running AIDO phenotype simulation."""
    st.info("Placeholder: Running AIDO phenotype simulation...")
    # In reality: Prepare input, run AIDO via subprocess, parse output.
    results = {
        'variant_id': filtered_variants_df['variant_id'].tolist(),
        'predicted_phenotype_effect': [f"Effect_{i}" for i in range(len(filtered_variants_df))]
    }
    st.success("Placeholder: AIDO simulation complete.")
    return pd.DataFrame(results)


# --- Streamlit App UI ---

st.set_page_config(page_title="Pulmonary Endothelium Toolkit", layout="wide")
st.title("🧬 Pulmonary Endothelium Toolkit")
st.markdown("""
Welcome to the Pulmonary Endothelium Toolkit (v2). This version includes functional
data fetching from uploaded VCF files (using `pysam`) and the GWAS Catalog API (using `requests`).
""")

# Initialize Database on first run
if 'db_initialized' not in st.session_state:
    init_db()
    st.session_state.db_initialized = True

# --- Sidebar Navigation ---
st.sidebar.header("Workflow Steps")
app_mode = st.sidebar.selectbox("Choose a section:",
    ["🏠 Home", "💾 Data Input & Processing", "📊 Database Query & Results", "ℹ️ About"]
)

# --- Main Page Content ---

if app_mode == "🏠 Home":
    st.header("Overview")
    st.markdown("""
    This toolkit aims to implement the computational aspects of the research program:
    * **Integrate Data:** Fetch variants from uploaded VCF files or the GWAS Catalog API.
    * **Analyze Variants:** Score variants (Evo2 - Placeholder), filter by cell type (ENCODE - Placeholder), simulate effects (AIDO - Placeholder).
    * **Store & Query:** Manage variants, scores, and annotations in a local database.
    * **Visualize:** Provide ways to explore the results.

    Use the sidebar to navigate through the different steps of the workflow.
    """)
    st.subheader("Project Focus Areas:")
    st.markdown("""
    - Endothelium in Lung Inflammation and Fibrosis (Oncology Population)
    - Hypoxia Regulation (GATA2-AS1, EPAS1)
    - Endothelial-Mesenchymal Transition (Radiation/Drug-induced)
    - Pulmonary Vascular Diseases (CTEPH, PAH)
    - Malignant Pleural Mesothelioma
    """)

elif app_mode == "💾 Data Input & Processing":
    st.header("Data Input and Pipeline Execution")

    # --- Step 1: Fetch Variants ---
    st.subheader("Step 1: Fetch or Upload Variants")
    data_source = st.radio("Select Data Source:", ["Upload VCF File", "Fetch from GWAS Catalog API"])

    variants_df = pd.DataFrame() # Initialize empty dataframe
    annotations_df = pd.DataFrame() # Initialize for GWAS annotations

    if data_source == "Upload VCF File":
        uploaded_file = st.file_uploader("Choose a VCF file (.vcf, .vcf.gz)", type=['vcf', 'gz'])
        if uploaded_file is not None:
            # Button to trigger parsing
            if st.button("Parse VCF File"):
                with st.spinner(f"Parsing {uploaded_file.name}..."):
                    variants_df = parse_vcf(uploaded_file)
                    st.session_state.variants_to_process = variants_df # Store in session state
                    st.session_state.annotations_to_process = pd.DataFrame() # No annotations from VCF parse

    elif data_source == "Fetch from GWAS Catalog API":
        trait = st.text_input("Enter GWAS Trait (e.g., pulmonary fibrosis, hypertension):", "pulmonary fibrosis")
        if st.button("Fetch GWAS Data"):
             with st.spinner(f"Fetching data for trait '{trait}' from GWAS Catalog..."):
                variants_df, annotations_df = fetch_gwas_catalog_by_trait(trait)
                st.session_state.variants_to_process = variants_df # Store variants
                st.session_state.annotations_to_process = annotations_df # Store annotations


    # --- Display fetched/uploaded variants ---
    if 'variants_to_process' in st.session_state and not st.session_state.variants_to_process.empty:
        st.write("Variants Loaded/Fetched:")
        st.dataframe(st.session_state.variants_to_process.head())

        if 'annotations_to_process' in st.session_state and not st.session_state.annotations_to_process.empty:
             st.write("Annotations Fetched (from GWAS):")
             st.dataframe(st.session_state.annotations_to_process.head())

        # --- Step 1b: Add Variants and Annotations to Database ---
        if st.button("Add Loaded Data to Database"):
            added_variants_count = 0
            skipped_variants_count = 0
            added_annotations_count = 0
            variant_db_ids_map = {} # Map variant_id string to DB id
            ids_for_scoring = []

            df_variants = st.session_state.variants_to_process
            df_annotations = st.session_state.annotations_to_process

            with st.spinner("Adding data to database..."):
                # Add Variants
                for index, row in df_variants.iterrows():
                    variant_id_str = row['variant_id']
                    # Use db_execute which handles UNIQUE constraint silently now
                    variant_db_id = db_execute("""
                        INSERT INTO variants (variant_id, chromosome, position, ref_allele, alt_allele, source_dataset, related_gene)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    """, (variant_id_str, row['chromosome'], row['position'], row['ref_allele'], row['alt_allele'], row['source_dataset'], row.get('related_gene')))

                    if variant_db_id:
                        added_variants_count += 1
                        variant_db_ids_map[variant_id_str] = variant_db_id
                        ids_for_scoring.append(variant_db_id)
                    else:
                        # If insert failed (likely duplicate), try to get existing ID
                        existing = db_query("SELECT id FROM variants WHERE variant_id = ?", (variant_id_str,))
                        if not existing.empty:
                            existing_id = existing['id'].iloc[0]
                            variant_db_ids_map[variant_id_str] = existing_id
                            ids_for_scoring.append(existing_id) # Add existing ID for potential re-scoring
                            skipped_variants_count += 1
                        else:
                             # Insert failed for other reason
                             st.warning(f"Failed to add or find variant: {variant_id_str}")


                # Add Annotations (if any) - requires mapping variant_id string to DB ID
                if not df_annotations.empty:
                    for index, row in df_annotations.iterrows():
                        variant_id_str = row['variant_id_ref'] # The linking ID from the fetch function
                        variant_db_id = variant_db_ids_map.get(variant_id_str) # Find the DB ID we stored/retrieved

                        if variant_db_id:
                            # Add check to avoid duplicate annotations for the same variant/type/source? Maybe later.
                            annot_id = db_execute("""
                                INSERT INTO annotations (variant_id, annotation_type, annotation_value, source, details)
                                VALUES (?, ?, ?, ?, ?)
                            """, (variant_db_id, row['annotation_type'], row['annotation_value'], row['source'], row.get('details')))
                            if annot_id:
                                added_annotations_count += 1
                        else:
                            st.warning(f"Could not find database ID for variant '{variant_id_str}' to add annotation.")

            st.success(f"Added {added_variants_count} new variants, skipped {skipped_variants_count} duplicates. Added {added_annotations_count} annotations.")
            # Remove duplicates before storing IDs for scoring
            st.session_state.variant_ids_for_scoring = list(set(ids_for_scoring))


    st.divider()

    # --- Step 2: Score Variants (Evo2 Placeholder) ---
    st.subheader("Step 2: Score Variants (Evo2 Placeholder)")
    if 'variant_ids_for_scoring' in st.session_state and st.session_state.variant_ids_for_scoring:
        st.write(f"Ready to score {len(st.session_state.variant_ids_for_scoring)} unique variants.")
        st.write(f"DB IDs: {st.session_state.variant_ids_for_scoring[:10]}...") # Show first 10

        # Retrieve variant details from DB needed for scoring input
        ids_tuple = tuple(st.session_state.variant_ids_for_scoring)
        # Ensure the tuple is not empty before querying
        if ids_tuple:
            variants_for_scoring_df = db_query(f"SELECT id, variant_id, chromosome, position FROM variants WHERE id IN ({','.join('?'*len(ids_tuple))})", ids_tuple)

            if not variants_for_scoring_df.empty:
                if st.button("Run Evo2 Scoring (Placeholder)"):
                    with st.spinner("Running Evo2 Scoring (Placeholder)..."):
                        scores_df = run_evo2_scoring(variants_for_scoring_df)
                        st.session_state.scores_to_process = scores_df # Store scores

                        # Add scores to database
                        added_scores = 0
                        if 'scores_to_process' in st.session_state and not st.session_state.scores_to_process.empty:
                            scores_df_to_add = st.session_state.scores_to_process
                            for index, row in scores_df_to_add.iterrows():
                                # Simple insert, could add checks for existing scores later
                                score_id = db_execute("""
                                    INSERT INTO scores (variant_id, tool, score_type, score_value, cell_type_context)
                                    VALUES (?, ?, ?, ?, ?)
                                """, (int(row['variant_id_internal']), row['tool'], row['score_type'], row['score_value'], row['cell_type_context']))
                                if score_id:
                                    added_scores += 1
                            st.success(f"Added {added_scores} scores to the database.")
                            st.dataframe(scores_df_to_add.head())
            else:
                st.warning("Could not retrieve variant details from database for scoring.")
        else:
            st.info("No variant IDs available for scoring.")

    else:
        st.info("Load/fetch variants and add them to the database first.")

    st.divider()

    # --- Step 3: Filter by Cell Type (ENCODE Placeholder) ---
    st.subheader("Step 3: Filter by Cell Type Context (ENCODE Placeholder)")
    # This part remains largely placeholder as ENCODE integration is complex
    # Check if scores exist in the database to be filtered
    scores_exist = db_query("SELECT COUNT(*) as count FROM scores")
    if not scores_exist.empty and scores_exist['count'].iloc[0] > 0:

        st.write("Scores available in database for filtering.")
        cell_type = st.selectbox("Select Target Cell Type:", ["ENDOS", "Fibroblast", "Epithelial", "Immune Cell"]) # Example cell types

        if st.button(f"Apply {cell_type} Filter (Placeholder)"):
            # Placeholder logic: Simulate filtering and potentially updating DB
            # In reality: Query scores, apply ENCODE logic, maybe update 'cell_type_context' in DB
            st.info("Simulating ENCODE filtering...")
            # Example: Query scores that might be relevant (e.g., context 'All')
            relevant_scores = db_query("SELECT id, variant_id FROM scores WHERE cell_type_context = 'All' LIMIT 10") # Example query
            if not relevant_scores.empty:
                 # Simulate updating some scores in the DB (replace with actual logic)
                 num_to_update = max(1, len(relevant_scores) // 2) # Update half
                 ids_to_update = tuple(relevant_scores['id'].head(num_to_update).tolist())
                 if ids_to_update:
                     update_count = db_execute(f"UPDATE scores SET cell_type_context = ? WHERE id IN ({','.join('?'*len(ids_to_update))})", (cell_type,) + ids_to_update)
                     # Note: db_execute returns lastrowid for INSERT, not useful for UPDATE count. Need separate query to confirm.
                     st.success(f"Placeholder: Simulated updating context to '{cell_type}' for up to {num_to_update} scores in DB.")
                     st.session_state.trigger_aido = True # Set flag for AIDO step
                 else:
                     st.warning("No scores found with context 'All' to simulate update.")

            else:
                 st.warning("No scores with context 'All' found in DB to simulate filtering.")

    else:
        st.info("Score variants first (Step 2).")

    st.divider()

    # --- Step 4: Simulate Phenotypes (AIDO Placeholder) ---
    st.subheader("Step 4: Simulate Phenotypic Effects (AIDO Placeholder)")
    # Trigger AIDO based on successful filtering simulation or manual button
    # Check if the flag was set in the previous step
    run_aido_flag = st.session_state.get('trigger_aido', False)

    if run_aido_flag:
        st.write("Filtered variants/scores potentially ready for simulation (based on placeholder filter).")
        # Query variants associated with the filtered scores (e.g., context = 'ENDOS')
        filtered_variants_for_aido = db_query("""
            SELECT v.id, v.variant_id
            FROM variants v JOIN scores s ON v.id = s.variant_id
            WHERE s.cell_type_context = ?
            LIMIT 50
        """, (cell_type,)) # Use cell_type from previous step selection

        if not filtered_variants_for_aido.empty:
             if st.button("Run AIDO Simulation (Placeholder)"):
                with st.spinner("Running AIDO Simulation (Placeholder)..."):
                    # Prepare input for AIDO (using variant_id strings)
                    aido_input_df = filtered_variants_for_aido[['variant_id']] # Pass df with 'variant_id' column
                    simulation_results = run_aido_simulation(aido_input_df)
                    st.success("Placeholder: AIDO simulation finished.")
                    st.dataframe(simulation_results.head())
                    # Reset flag
                    st.session_state.trigger_aido = False
        else:
             st.warning(f"No variants found with score context '{cell_type}' in DB for AIDO simulation.")
             # Reset flag
             st.session_state.trigger_aido = False

    else:
        st.info("Apply cell type filter first (Step 3) to enable AIDO simulation.")


elif app_mode == "📊 Database Query & Results":
    # --- Database Query Section (Largely unchanged but uses updated schema if needed) ---
    st.header("Database Query and Results")

    st.subheader("Query Variants")
    # ... (query interface remains the same as v1) ...
    query_gene = st.selectbox("Filter by Gene (if annotated):", ["All"] + GENES_OF_INTEREST)
    query_source = st.selectbox("Filter by Source Dataset:", ["All", "VCF Upload", "GWAS Catalog API"]) # Updated sources
    query_variant_id = st.text_input("Filter by Variant ID (e.g., rs123, chr1:12345):")

    sql = "SELECT v.id, v.variant_id, v.chromosome, v.position, v.ref_allele, v.alt_allele, v.source_dataset, v.related_gene, v.date_added FROM variants v WHERE 1=1"
    params = []

    if query_gene != "All":
        st.warning("Filtering by gene requires 'related_gene' to be populated (e.g., via external annotation step).")
        sql += " AND v.related_gene = ?"
        params.append(query_gene)
    if query_source != "All":
        # Handle potential naming variations if needed
        source_like = f"%{query_source}%" if query_source == "VCF Upload" else query_source
        sql += " AND v.source_dataset LIKE ?" if query_source == "VCF Upload" else " AND v.source_dataset = ?"
        params.append(source_like)
    if query_variant_id:
        sql += " AND v.variant_id LIKE ?"
        params.append(f"%{query_variant_id}%")

    sql += " ORDER BY v.date_added DESC LIMIT 100"

    variants_results = db_query(sql, tuple(params))
    st.write(f"Found {len(variants_results)} variants matching criteria:")
    st.dataframe(variants_results)

    st.divider()

    st.subheader("Query Scores")
    # ... (query interface remains the same as v1) ...
    score_min = st.slider("Minimum Score Value:", 0.0, 1.0, 0.0, 0.01)
    score_max = st.slider("Maximum Score Value:", 0.0, 1.0, 1.0, 0.01)
    score_tool = st.selectbox("Filter by Scoring Tool:", ["All", "Evo2_Placeholder"])
    score_cell_type = st.selectbox("Filter by Cell Type Context:", ["All", "ENDOS", "Fibroblast", "Epithelial", "Immune Cell", "Other"]) # Added more options

    sql_scores = """
        SELECT s.id, v.variant_id, s.tool, s.score_type, s.score_value, s.cell_type_context, s.date_scored
        FROM scores s
        JOIN variants v ON s.variant_id = v.id
        WHERE s.score_value BETWEEN ? AND ?
    """
    params_scores = [score_min, score_max]

    if score_tool != "All":
        sql_scores += " AND s.tool = ?"
        params_scores.append(score_tool)
    if score_cell_type != "All":
        sql_scores += " AND s.cell_type_context = ?"
        params_scores.append(score_cell_type)

    sql_scores += " ORDER BY s.score_value DESC LIMIT 100"

    scores_results = db_query(sql_scores, tuple(params_scores))
    st.write(f"Found {len(scores_results)} scores matching criteria:")
    st.dataframe(scores_results)


    st.divider()

    st.subheader("Query Annotations")
    # ... (query interface remains the same as v1, but reflects GWAS data) ...
    annot_type = st.selectbox("Filter by Annotation Type:", ["All", "Reported Trait", "Disease Association", "Pathway", "Gene Function"]) # Added Reported Trait
    annot_value = st.text_input("Filter by Annotation Value (e.g., pulmonary fibrosis, CTEPH):")
    annot_source = st.selectbox("Filter by Annotation Source:", ["All", "GWAS Catalog API", "Manual Curation"])

    sql_annots = """
        SELECT a.id, v.variant_id, a.annotation_type, a.annotation_value, a.source, a.details, a.date_annotated
        FROM annotations a
        JOIN variants v ON a.variant_id = v.id
        WHERE 1=1
    """
    params_annots = []

    if annot_type != "All":
        sql_annots += " AND a.annotation_type = ?"
        params_annots.append(annot_type)
    if annot_value:
        sql_annots += " AND a.annotation_value LIKE ?"
        params_annots.append(f"%{annot_value}%")
    if annot_source != "All":
        sql_annots += " AND a.source LIKE ?" # Use LIKE for flexibility (e.g., matching 'GWAS Catalog API (ACCESSION)')
        params_annots.append(f"%{annot_source}%")


    sql_annots += " ORDER BY a.date_annotated DESC LIMIT 100"

    annots_results = db_query(sql_annots, tuple(params_annots))
    st.write(f"Found {len(annots_results)} annotations matching criteria:")
    st.dataframe(annots_results)


elif app_mode == "ℹ️ About":
    st.header("About This Toolkit (v2)")
    st.markdown("""
    This Streamlit application provides a framework for the computational workflow outlined in the
    "Decoding Endothelial Influence on Pulmonary Disease" research program.

    **Version 2 Updates:**
    - **Functional Data Input:** Replaced placeholders with functional code:
        - Parses uploaded VCF files using `pysam`.
        - Fetches data from the GWAS Catalog API using `requests`.
    - **Database Population:** Adds fetched/parsed variants and annotations to the SQLite database.
    - **Improved DB Handling:** Silently handles duplicate variant entries during insertion.

    **Remaining Placeholders:**
    - Evo2 scoring (`run_evo2_scoring`).
    - ENCODE filtering (`filter_encode_data`).
    - AIDO simulation (`run_aido_simulation`).
    - Gene annotation based on variant position/ID.

    **Dependencies Added:**
    - `pysam`: For VCF parsing (`pip install pysam`).
    - `requests`: For API calls (`pip install requests`).

    **Future Development:**
    - Integrate actual command-line execution of Evo2, AIDO, etc. using `subprocess`.
    - Implement ENCODE data fetching/filtering logic.
    - Add gene annotation step (e.g., using external tools/APIs like VEP or querying Ensembl/NCBI).
    - Develop visualization modules (Plotly).
    - Package the application properly.
    """)
    st.info(f"Database location: {os.path.abspath(DB_PATH)}")

