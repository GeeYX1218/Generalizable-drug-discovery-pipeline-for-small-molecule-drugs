```markdown
# 🧬 Integrated Computational Drug Discovery Pipeline
### Generalized for Small Molecule Inhibitors

## 📌 Overview
This Jupyter Notebook implements an end-to-end *in silico* drug discovery pipeline. It automates the process of retrieving bioactivity data, analyzing Structure-Activity Relationships (SAR), generating novel drug candidates, and validating them via molecular docking, machine learning and Multiparameter optimization.

**🚀 Pre-execution Configuration Checklist**
**1. Target Definition (Step 1 - Mandatory)**
[ ] Define Target ChEMBL ID: Have you found the Single Protein ID (e.g., CHEMBL2842 for mTOR)?
    Action: Set TARGET_CHEMBL_ID. (Recommended over searching by name to avoid "dirty" data).

[ ] Set Search Term: Does this match your target protein's name?
    Action: Set TARGET_SEARCH_TERM (e.g., 'EGFR', 'KRAS').

[ ] Set Project Name: Is this short and filesystem-friendly (no spaces/special chars)?
    Action: Set TARGET_NAME (e.g., 'egfr', 'braf').

[ ] Configure Mutant Filter: Are you targeting a specific mutation (e.g., G12D, T790M) or Wild Type?
    Action: Set MUTANT_FILTER to a string or None.

[ ] Set Potency Threshold: Is 50.0 nM appropriate, or is your target "tough" (requiring 100.0 or 1000.0 nM)?
    Action: Adjust POTENCY_CUTOFF_NM.

**2. External Access (Step 3 - Recommended)**
[ ] Set NCBI Email: Have you entered a valid email address to prevent NCBI from blocking your literature search?
    Action: Update Entrez.email.

[ ] Select PDB Strategy: Will you rely on Auto-Detection (riskier) or Manual Override (safer)?
    Action: If Manual, set MANUAL_PDB_ID (e.g., '4JT6').

**3. Docking & Structure (Step 5 - Critical)**
[ ] ⚠️ The "Belly Button" Safety Check: Does your chosen PDB have a known co-crystallized ligand?
    If YES: You are safe. The pipeline will center the grid on the ligand.
    If NO (or Unsure): You must find the Residue Number of the active site.
    Action: Have this number ready to enter into TARGET_RESIDUE_ID in the else block (~line 130) if the "Using geometric center" warning appears.

[ ] Tuning Docking Parameters: Do you need a larger search space or higher precision?
    Action: Review BOX_SIZE (Default: 20.0) and EXHAUSTIVENESS (Default: 8).

**✅ Ready to Run?**
If you have checked all the boxes above, you can confidently run the cells. Remember to watch the output of Step 5 closely for the "Geometric Center" warning!**

## ⚙️ Configuration Guide (Crucial)
**To switch targets (e.g., to EGFR or BRAF), you must modify specific lines in the code.** Below is the line-by-line analysis of where user input is required.

### 🔹 Step 1: Data Acquisition (Mandatory)
This is the most critical setup block. Changing variables here propagates through the entire pipeline.
*Located in the second code cell, lines ~25-45.*

```python
# ==============================================================================
# 🛠️ USER CONFIGURATION (CHANGE THIS BLOCK TO SWITCH TARGETS)
# ==============================================================================

# 1. Target ID (Strongly Recommended to include): If you know the exact ChEMBL ID (e.g., CHEMBL203 for EGFR).
#    Leave as None to search by name.
TARGET_CHEMBL_ID = None or "CHEMBL203" 

# 2. Search Term: The protein name to query in ChEMBL.
#    CHANGE THIS: Replace 'KRAS' with 'EGFR', 'BRAF', etc.
TARGET_SEARCH_TERM = 'KRAS'

# 3. Project Name: Used for naming all output files (e.g., 'egfr_inhibitors.csv').
#    CHANGE THIS: Replace 'kras' with 'egfr', 'her2', etc.
TARGET_NAME = 'kras'

# 4. Mutant Filter: Filters data for specific mutations (e.g., 'T790M' for EGFR).
#    CHANGE THIS: Set to None for Wild Type, or specific mutation code.
MUTANT_FILTER = "G12D" or None

# 5. Potency Threshold: Only keeps inhibitors with IC50 better than this value.
#    ADJUST IF NEEDED: 50.0 nM is standard for "potent" hits.
POTENCY_CUTOFF_NM = 50.0 

```

### 🔹 Step 3: Genomic & Structural Analysis (Optional)

```python
# PDB Structure: The pipeline tries to auto-detect a PDB by searching with ChEMBL and UniProt. 
# IF AUTO-DETECT FAILS or you want a specific crystal structure:
# 🛠️ MANUAL PDB OVERRIDE
# Set this to a string (e.g., "4JT6") to force a specific structure.
# Set to None to enable automatic searching.
# ------------------------------------------------------------------------------
MANUAL_PDB_ID = "7RT1" # default None (Automatic search), unless the searching result doesn't match mutant type/required structure

```
```python
# Email for NCBI Entrez: Required for PubMed literature mining.
# RECOMMENDATION: Change to your actual email to avoid NCBI rate limits/blocks.
Entrez.email = "your.email@nus.edu.sg" 
```

### 🔹 Step 5: Molecular Docking (Virtual Screening) (Optional)
# Docking Parameters:
# ADJUST IF NEEDED: Increase for higher precision, decrease for speed.
BOX_SIZE = 20.0       # Size of the search space (Angstroms)
EXHAUSTIVENESS = 8    # How hard Vina searches for poses (1-8 is standard)

```

# ⚠️ Critical Warning (The "Belly Button" Problem): If your chosen PDB does not have a co-crystallized ligand, the pipeline defaults to the Geometric Center, which is often inaccurate (targeting the protein core).

# Action Required: If you see the warning "Using geometric center" in the output, go to the fallback else block (~line 130) and manually set TARGET_RESIDUE_ID to the residue number of your active site (e.g., 12 for KRAS).

```python
 # Calculate binding site center from ligand coordinates
    if ligand_coords:
        ligand_coords = np.array(ligand_coords)
        center = ligand_coords.mean(axis=0)
        center_x, center_y, center_z = center
        print(f"✓ Co-crystallized ligand found: {len(ligand_coords)} atoms")
        print(f"✓ Binding site center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
    else:
        # BETTER FALLBACK: Target a specific Active Site Residue
        # Change '12' to the residue number of your active site (e.g., G12 in KRAS)
        TARGET_RESIDUE_ID = 12  
        
        print(f"⚠ Warning: Ligand not found. Targeting Residue {TARGET_RESIDUE_ID}...")
```
---

## 🚀 Pipeline Workflow & Example Outputs

### **Step 0: Dependencies & Environment**

* **Action:** Installs `rdkit`, `openbabel`, `autodock-vina`, and `chembl_webresource_client`.
* **Note:** Must be run once at the start of the session.

### **Step 1: Data Acquisition**

Retrieves experimental bioactivity data from ChEMBL, filters for high potency (), cleans chemical structures, and removes duplicates.

* **Outputs:**
* `kras_inhibitors_cleaned.csv`: The curated dataset of potent inhibitors.
* `kras_pic50_distribution.png`: Histogram showing the potency spread of the dataset.



### **Step 2: SAR & Scaffold Analysis**

Performs Bemis-Murcko scaffold decomposition to identify "privileged scaffolds" (core structures that yield high potency).

* **Outputs:**
* `kras_scaffold_analysis.csv`: List of scaffolds ranked by frequency and potency.
* `kras_sar_analysis.csv`: Extended dataset with molecular descriptors (MW, LogP, TPSA).
* `kras_chemical_space_pca.png`: PCA plot visualizing chemical diversity.
* `kras_physicochemical_properties.png`: Violin plots of Lipinski properties.
* `kras_top_scaffolds_potency.png`: Bar chart of the best-performing scaffolds.



### **Step 3: Genomic & Structural Context**

Searches RCSB PDB for available crystal structures and PubMed for literature on resistance mechanisms.

* **Outputs:**
* `kras_structural_analysis.json`: Metadata on the best available PDB structure (resolution, ligands).
* `kras_literature_findings.txt`: Abstracts of relevant papers (e.g., "KRAS G12D inhibitors").



### **Step 4: Generative Design**

Uses biochemical mutation strategies (bioisosteric replacement, ring modification) on the top scaffolds from Step 2 to generate novel analogs. Filters for drug-likeness (QED) and Synthetic Accessibility (SA).

* **Outputs:**
* `kras_selected_seeds.csv`: The parent molecules used for generation.
* `kras_generated_candidates.csv`: All valid generated molecules.
* `kras_top20_generated_candidates.csv`: The top 20 candidates ranked by combined score.
* `kras_top20_structures.png`: Grid image of the top 20 novel structures.
* `kras_generation_pca/tsne.png`: Plots showing where new candidates sit in chemical space relative to known inhibitors.



### **Step 5: Virtual Screening (Docking)**

Prepares the receptor (PDB) and ligands (PDBQT), then runs AutoDock Vina to predict binding affinity (kcal/mol).

* **Outputs:**
* `kras_docking_results.csv`: Ranked list of candidates by binding affinity.
* `kras_docking_scores.png`: Bar chart of binding energies.



### **Step 6: Machine Learning Prediction**

Trains a Random Forest regressor on the ChEMBL data (from Step 1) to predict the  of the generated candidates (from Step 4) as a secondary validation.

* **Outputs:**
* `kras_candidate_predictions.csv`: ML-predicted potency values.
* `model_performance.png`: Actual vs. Predicted plot for the test set.
* `ml_docking_comparison.png`: Correlation plot between ML predictions and Docking scores.



### **Step 7: ADME/Tox & Final Selection**

Evaluates pharmacokinetic properties (Absorption, Distribution, Metabolism, Excretion) and toxicity risks.

* **Outputs:**
* `kras_admet_analysis.csv`: Predicted ADMET properties.
* `kras_final_candidates.csv`: The final "Hit" list after all filters.
* `kras_candidate_radar_plot.png`: Multi-parameter visualization of the top candidates.



---

## ⚠️ Troubleshooting & Notes

1. **Standalone Execution:** If you run steps out of order (e.g., Step 4 without running Step 2), the script attempts to load files from the `results/` folder. Ensure the filename prefixes match your `TARGET_NAME`.
2. **PDB Download:** If Step 5 fails to download a PDB, manually upload a `.pdb` file to the `workflow/data` folder and the script will auto-detect it.
3. **Dependencies:** If you see `ModuleNotFoundError`, re-run **Step 0**.

```

```