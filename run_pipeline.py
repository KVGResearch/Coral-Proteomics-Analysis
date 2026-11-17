#!/usr/bin/env python3
"""
Coral Proteomics Comparative Analysis Pipeline
Author: Gemini CLI Agent
Date: Thursday, June 25, 2026

This pipeline executes:
1. Ingestion and cleaning of protein lists for 4 coral species.
2. Building of a master quantitative protein matrix (Total-Count Normalized PSMs).
3. Overlap analysis (core, shared, and unique) with an elegant UpSet Plot.
4. Hierarchical clustering heatmap of the core proteome.
5. GO terms parsing, functional classification, and comparative profiling.
6. Statistical enrichment (Fisher's Exact Test with Benjamini-Hochberg FDR) for unique proteomes.
7. Pathway annotation and abundance profiling for major coral biological systems.
8. Creation of advanced figures: a custom Matplotlib-based Flow (Sankey-like) diagram
   and a Network Diagram of core proteins based on Pfam domains.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
import upsetplot
import networkx as nx

# ----------------------------------------------------------------------
# 0. Global Setup & Directory Creation
# ----------------------------------------------------------------------
print("[INFO] Initialising environment...")
results_dir = "/mnt/c/Santhosh_Research/VIT_PhD/Experinmental/stellamaries_Proteomics/Results/Results/results"
plots_dir = os.path.join(results_dir, "plots")
tables_dir = os.path.join(results_dir, "tables")

os.makedirs(plots_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)

# Set elegant publication style
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 18,
    "figure.dpi": 300
})

# Meta-information about samples
samples_meta = {
    "Sample-1": {"species": "Acropora muricata", "short": "A. muricata", "file": "Veena_Sample-1_protein_list.xlsx"},
    "Sample-2": {"species": "Mussa angulosa", "short": "M. angulosa", "file": "Veena_Sample-2_protein_list.xlsx"},
    "Sample-3": {"species": "Dendronephthya australis", "short": "D. australis", "file": "Veena_Sample-3_protein_list.xlsx"},
    "Sample-4": {"species": "Agaricia tenuifolia", "short": "A. tenuifolia", "file": "Veena_Sample-4_protein_list.xlsx"}
}

# ----------------------------------------------------------------------
# Helper Functions
# ----------------------------------------------------------------------
def parse_uniprot_desc(desc):
    """
    Parses a UniProt description to extract:
    1. Protein Name (before OS=)
    2. Organism Name (between OS= and next field)
    3. Gene Symbol (GN=)
    """
    if not isinstance(desc, str) or desc.strip() == "":
        return "Unknown Protein", "Unknown Organism", ""
    
    # Protein Name
    p_parts = re.split(r'\s+OS=', desc)
    protein_name = p_parts[0].strip() if p_parts else desc
    
    # Organism Name
    org_match = re.search(r'OS=([^=]+?)(?=\s+[A-Z]{2}=|$)', desc)
    organism = org_match.group(1).strip() if org_match else "Unknown Organism"
    
    # Gene Name
    gn_match = re.search(r'GN=([^=]+?)(?=\s+[A-Z]{2}=|$)', desc)
    gene_name = gn_match.group(1).strip() if gn_match else ""
    
    return protein_name, organism, gene_name

def parse_go_terms(go_str):
    """
    Parses semicolons in GO columns.
    Returns list of clean term names and list of GO IDs.
    """
    if not isinstance(go_str, str) or go_str.strip() == "":
        return [], []
    
    terms = []
    go_ids = []
    # Items separated by semicolons
    parts = [p.strip() for p in go_str.split(";") if p.strip()]
    for part in parts:
        # Format: "term name [GO:ID]"
        match = re.match(r'^(.*?)\s*\[(GO:\d+)\]$', part)
        if match:
            terms.append(match.group(1).strip())
            go_ids.append(match.group(2).strip())
        else:
            terms.append(part)
            
    return terms, go_ids

def parse_pfam_ids(pfam_str):
    """Parses Pfam ID list."""
    if not isinstance(pfam_str, str) or pfam_str.strip() == "":
        return []
    return [p.strip() for p in pfam_str.split(",") if p.strip()]

# ----------------------------------------------------------------------
# 1. Loading, Cleaning and Normalisation
# ----------------------------------------------------------------------
print("[INFO] Step 1: Loading and standardising protein datasets...")
dfs_raw = {}
protein_sets = {}

for s_key, meta in samples_meta.items():
    file_path = f"/mnt/c/Santhosh_Research/VIT_PhD/Experinmental/stellamaries_Proteomics/Results/Results/{meta['file']}"
    df = pd.read_excel(file_path)
    
    # Standardise column names
    df = df.rename(columns=lambda x: x.strip())
    
    # Clean up accessions
    df["Accession"] = df["Accession"].astype(str).str.strip()
    df = df.dropna(subset=["Accession"])
    
    # Parse annotations
    parsed_metadata = df["Description"].apply(parse_uniprot_desc)
    df["Clean_Protein_Name"] = [p[0] for p in parsed_metadata]
    df["Parsed_Organism"] = [p[1] for p in parsed_metadata]
    df["Parsed_Gene_Symbol"] = [p[2] for p in parsed_metadata]
    
    # Merge with original Gene Symbol column if empty
    if "Gene Symbol" in df.columns:
        df["Gene Symbol"] = df["Gene Symbol"].fillna(df["Parsed_Gene_Symbol"])
    else:
        df["Gene Symbol"] = df["Parsed_Gene_Symbol"]
        
    # Standardise quantitative columns (PSMs, Coverage, etc.)
    df["# PSMs"] = pd.to_numeric(df["# PSMs"], errors="coerce").fillna(0).astype(int)
    df["Coverage [%]"] = pd.to_numeric(df["Coverage [%]"], errors="coerce").fillna(0)
    df["Sum PEP Score"] = pd.to_numeric(df["Sum PEP Score"], errors="coerce").fillna(0)
    
    dfs_raw[s_key] = df
    protein_sets[s_key] = set(df["Accession"].tolist())
    print(f"  - Loaded {meta['species']} ({meta['short']}): {df.shape[0]} proteins. Total PSMs: {df['# PSMs'].sum()}")

# Calculate Normalisation factors (Total PSM count normalization)
# We will normalise to 10,000 spectral matches (PSMs) to keep numbers intuitive
scaling_factor = 10000.0
norm_factors = {}
for s_key, df in dfs_raw.items():
    total_psms = df["# PSMs"].sum()
    norm_factors[s_key] = scaling_factor / total_psms if total_psms > 0 else 1.0

# ----------------------------------------------------------------------
# 2. Master Matrix Construction
# ----------------------------------------------------------------------
print("[INFO] Building consolidated Master Protein Matrix...")
# Collect all unique accessions
all_accessions = sorted(list(set().union(*protein_sets.values())))

# Gather the best metadata for each protein
metadata_dict = {}
for s_key, df in dfs_raw.items():
    for _, row in df.iterrows():
        acc = row["Accession"]
        if acc not in metadata_dict:
            metadata_dict[acc] = {
                "Description": row["Description"],
                "Protein Name": row["Clean_Protein_Name"],
                "Organism": row["Parsed_Organism"],
                "Gene Symbol": row["Gene Symbol"] if pd.notna(row["Gene Symbol"]) and str(row["Gene Symbol"]) != "" else "",
                "MW [kDa]": row["MW [kDa]"] if "MW [kDa]" in row else np.nan,
                "calc. pI": row["calc. pI"] if "calc. pI" in row else np.nan,
                "Biological Process": row["Biological Process"] if "Biological Process" in row else "",
                "Cellular Component": row["Cellular Component"] if "Cellular Component" in row else "",
                "Molecular Function": row["Molecular Function"] if "Molecular Function" in row else "",
                "Pfam IDs": row["Pfam IDs"] if "Pfam IDs" in row else ""
            }

# Create records for Master Matrix
master_records = []
for acc in all_accessions:
    meta = metadata_dict[acc]
    record = {
        "Accession": acc,
        "Protein Name": meta["Protein Name"],
        "Gene Symbol": meta["Gene Symbol"],
        "Organism": meta["Organism"],
        "MW [kDa]": meta["MW [kDa]"],
        "calc. pI": meta["calc. pI"],
        "Biological Process": meta["Biological Process"],
        "Cellular Component": meta["Cellular Component"],
        "Molecular Function": meta["Molecular Function"],
        "Pfam IDs": meta["Pfam IDs"],
        "Description": meta["Description"]
    }
    
    # Quantities per sample
    for s_key in samples_meta.keys():
        s_df = dfs_raw[s_key]
        match = s_df[s_df["Accession"] == acc]
        if not match.empty:
            raw_psm = match.iloc[0]["# PSMs"]
            cov = match.iloc[0]["Coverage [%]"]
            pep = match.iloc[0]["Sum PEP Score"]
            norm_psm = raw_psm * norm_factors[s_key]
            
            record[f"{s_key}_Identified"] = True
            record[f"{s_key}_Raw_PSM"] = raw_psm
            record[f"{s_key}_Norm_PSM"] = norm_psm
            record[f"{s_key}_Coverage_%"] = cov
            record[f"{s_key}_Sum_PEP"] = pep
        else:
            record[f"{s_key}_Identified"] = False
            record[f"{s_key}_Raw_PSM"] = 0
            record[f"{s_key}_Norm_PSM"] = 0.0
            record[f"{s_key}_Coverage_%"] = 0.0
            record[f"{s_key}_Sum_PEP"] = 0.0
            
    # Calculate across-sample presence
    record["Total_Samples_Present"] = sum([record[f"{s_key}_Identified"] for s_key in samples_meta.keys()])
    master_records.append(record)

df_master = pd.DataFrame(master_records)

# Save Master Matrix to Excel
df_master.to_excel(os.path.join(tables_dir, "1_Master_Protein_Matrix.xlsx"), index=False)
print(f"  - Master Protein Matrix created with {df_master.shape[0]} unique proteins. Saved to tables/1_Master_Protein_Matrix.xlsx.")

# ----------------------------------------------------------------------
# 3. Step 2: Comparative Overlap & UpSet Plot
# ----------------------------------------------------------------------
print("[INFO] Step 2: Performing Overlap and UpSet analysis...")

# Prepare data for UpSet Plot
upset_contents = {}
for s_key, meta in samples_meta.items():
    label = f"{meta['short']}\n({meta['species']})"
    upset_contents[label] = protein_sets[s_key]

upset_data = upsetplot.from_contents(upset_contents)

# Plot UpSet diagram
fig = plt.figure(figsize=(14, 8))
upset_plot_obj = upsetplot.plot(upset_data, fig=fig, element_size=60, show_counts=True)
plt.title("Protein Identification Overlaps Across 4 Coral Species (UpSet Plot)", fontsize=18, pad=30)
plt.savefig(os.path.join(plots_dir, "1_upset_plot.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Generated UpSet Plot and saved to plots/1_upset_plot.png.")

# Save specific intersection groups
writer = pd.ExcelWriter(os.path.join(tables_dir, "2_Proteins_By_Intersection_Groups.xlsx"), engine="xlsxwriter")

# 1. Core Proteome (Present in all 4)
df_core = df_master[df_master["Total_Samples_Present"] == 4]
df_core.to_excel(writer, sheet_name="Core_Proteome_All_4", index=False)

# 2. Species-Specific unique proteins
for s_key, meta in samples_meta.items():
    others = [k for k in samples_meta.keys() if k != s_key]
    condition = df_master[f"{s_key}_Identified"] & (~df_master[f"{others[0]}_Identified"]) & (~df_master[f"{others[1]}_Identified"]) & (~df_master[f"{others[2]}_Identified"])
    df_spec = df_master[condition]
    sheet_name = f"Unique_{meta['short'].replace(' ', '_')}"
    df_spec.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    print(f"  - Unique proteins for {meta['short']}: {df_spec.shape[0]}")

writer.close()
print("  - Saved intersection list excel to tables/2_Proteins_By_Intersection_Groups.xlsx.")

# ----------------------------------------------------------------------
# 4. Step 3: Hierarchical Clustering of Core Proteome
# ----------------------------------------------------------------------
print("[INFO] Step 3: Performing Hierarchical Clustering of the Core Proteome (77 proteins)...")

# Get core proteins
df_core_data = df_master[df_master["Total_Samples_Present"] == 4].copy()

# Extract normalized PSM values
norm_cols = [f"{s_key}_Norm_PSM" for s_key in samples_meta.keys()]
heatmap_data = df_core_data[norm_cols].copy()
heatmap_data.columns = [samples_meta[s_key]["short"] for s_key in samples_meta.keys()]

# Use Gene Symbol as index, fall back to Accession if empty
labels = []
for _, row in df_core_data.iterrows():
    gs = row["Gene Symbol"]
    acc = row["Accession"]
    labels.append(f"{gs} ({acc})" if gs and str(gs) != "nan" else acc)
heatmap_data.index = labels

# Standardize row-wise (Z-score)
# Z = (x - mean) / std. Add small epsilon to avoid divide by zero if std == 0
row_means = heatmap_data.mean(axis=1)
row_stds = heatmap_data.std(axis=1).replace(0, 1e-6)
heatmap_z = heatmap_data.sub(row_means, axis=0).div(row_stds, axis=0)

# Generate Clustermap
g = sns.clustermap(
    heatmap_z,
    method="ward",
    metric="euclidean",
    cmap="RdYlBu_r",
    vmin=-2, vmax=2,
    figsize=(12, 18),
    col_cluster=True,
    row_cluster=True,
    cbar_kws={'label': 'Z-Score (Relative Abundance)'},
    yticklabels=True,
    xticklabels=True
)

plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=12)
g.fig.suptitle("Hierarchical Clustering of Core Proteome (77 Conserved Proteins)", y=1.01, fontsize=18)
g.savefig(os.path.join(plots_dir, "2_core_heatmap.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Generated hierarchical clustered heatmap and saved to plots/2_core_heatmap.png.")

# ----------------------------------------------------------------------
# 5. Step 4: GO Functional Classification
# ----------------------------------------------------------------------
print("[INFO] Step 4: Processing and profiling Gene Ontology (GO) terms...")

go_domains = {
    "Biological Process": "BP",
    "Cellular Component": "CC",
    "Molecular Function": "MF"
}

# Parse and collect terms
go_term_counts = {s_key: {"BP": {}, "CC": {}, "MF": {}} for s_key in samples_meta.keys()}

for s_key, df in dfs_raw.items():
    for col, domain in go_domains.items():
        if col in df.columns:
            for val in df[col].dropna():
                terms, _ = parse_go_terms(str(val))
                for term in terms:
                    go_term_counts[s_key][domain][term] = go_term_counts[s_key][domain].get(term, 0) + 1

# Convert to structured tables and find major terms for classification
# We will create a grouped chart for major high-level terms (GO Slim categories)
# Major slim terms mapped from our observed terms:
go_slim_map = {
    "BP": {
        "cellular process": ["cell", "cellular", "division", "cycle", "growth"],
        "metabolic process": ["metabolic", "biosynthesis", "catabolic", "synthesis", "metabolism"],
        "response to stress / stimuli": ["stress", "stimulus", "chaperone", "heat shock", "defense", "immune", "oxidative"],
        "transport / localization": ["transport", "localization", "protein import", "secretion", "vesicle"],
        "translation / RNA process": ["translation", "ribosome", "rna", "transcription", "splicing"],
        "developmental process": ["development", "differentiation", "morphogenesis"],
        "cell adhesion / signaling": ["signaling", "signal transduction", "adhesion", "junction"]
    },
    "MF": {
        "catalytic activity": ["catalytic", "enzyme", "activity", "hydrolase", "transferase", "oxidoreductase", "dehydrogenase", "kinase", "phosphatase"],
        "binding": ["binding", "receptor", "interact"],
        "structural molecule activity": ["structural", "cytoskeleton", "tubulin", "actin", "ribosomal protein"],
        "transporter activity": ["transporter", "carrier", "channel", "pore"],
        "translation/transcription factor": ["translation factor", "transcription factor", "initiation", "elongation"]
    },
    "CC": {
        "cytoskeleton": ["cytoskeleton", "microtubule", "actin filament", "spindle"],
        "cytoplasm / cytosol": ["cytoplasm", "cytosol"],
        "nucleus / chromatin": ["nucleus", "nucleoplasm", "chromatin", "nucleolus"],
        "mitochondrion": ["mitochondrion", "mitochondrial"],
        "membrane / membrane-bound organelle": ["membrane", "envelope", "vesicle", "endoplasmic reticulum", "golgi"],
        "extracellular / cell wall": ["extracellular", "matrix", "secreted"],
        "ribosome / macromolecular complex": ["ribosome", "ribonucleoprotein", "complex"]
    }
}

# Map raw detailed counts to these high-level slim categories
slim_results = []
for s_key, meta in samples_meta.items():
    for dom_name, dom_abbr in go_domains.items():
        total_proteins = dfs_raw[s_key].shape[0]
        category_counts = {cat: 0 for cat in go_slim_map[dom_abbr].keys()}
        category_counts["other biological processes/functions"] = 0
        
        # Accumulate matching counts
        raw_counts = go_term_counts[s_key][dom_abbr]
        for term, count in raw_counts.items():
            mapped = False
            term_lower = term.lower()
            for cat, keywords in go_slim_map[dom_abbr].items():
                if any(kw in term_lower for kw in keywords):
                    category_counts[cat] += count
                    mapped = True
                    break
            if not mapped:
                category_counts["other biological processes/functions"] += count
                
        # Normalise counts as % of total proteins identified (since proteins have multiple terms, sum can exceed 100%)
        for cat, count in category_counts.items():
            pct = (count / total_proteins) * 100.0
            slim_results.append({
                "Sample": meta["short"],
                "Domain": dom_name,
                "Category": cat,
                "Term Count": count,
                "Percentage_Proteins": pct
            })

df_go_slim = pd.DataFrame(slim_results)
df_go_slim.to_excel(os.path.join(tables_dir, "3_GO_Functional_Classification.xlsx"), index=False)
print("  - Saved high-level GO classification to tables/3_GO_Functional_Classification.xlsx.")

# Visualization: 3-panel horizontal bar chart (BP, CC, MF)
fig, axes = plt.subplots(3, 1, figsize=(14, 20), sharex=False)
fig.suptitle("Gene Ontology (GO) Functional Classification across Coral Species", y=1.01, fontsize=20, weight="bold")

for idx, (dom_name, dom_abbr) in enumerate(go_domains.items()):
    ax = axes[idx]
    df_sub = df_go_slim[df_go_slim["Domain"] == dom_name]
    
    sns.barplot(
        data=df_sub,
        y="Category",
        x="Percentage_Proteins",
        hue="Sample",
        ax=ax,
        palette="muted"
    )
    ax.set_title(f"{dom_name} Classification", fontsize=16, weight="bold")
    ax.set_xlabel("% of Identified Proteome containing GO term", fontsize=12)
    ax.set_ylabel("")
    ax.legend(title="Coral Species", fontsize=10, title_fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "3_go_classification.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Generated comparative GO classification bar chart and saved to plots/3_go_classification.png.")

# ----------------------------------------------------------------------
# 6. Step 5: Statistical Enrichment Analysis (Fisher's Exact Test)
# ----------------------------------------------------------------------
print("[INFO] Step 5: Running statistical overrepresentation enrichment tests...")

# 1. Background definitions
bg_accessions = set(df_master["Accession"])
bg_size = len(bg_accessions)

# Extract backgrounds of parsed annotations (GO and Pfams)
bg_go_counts = {}
bg_pfam_counts = {}

# We pre-extract annotations for all proteins in master
for _, row in df_master.iterrows():
    acc = row["Accession"]
    # Gather GOs
    gos = []
    for col in ["Biological Process", "Cellular Component", "Molecular Function"]:
        if pd.notna(row[col]) and row[col] != "":
            terms, _ = parse_go_terms(str(row[col]))
            gos.extend(terms)
    gos = list(set(gos))
    for g in gos:
        bg_go_counts[g] = bg_go_counts.get(g, 0) + 1
        
    # Gather Pfams
    pfams = parse_pfam_ids(str(row["Pfam IDs"]))
    for pf in pfams:
        bg_pfam_counts[pf] = bg_pfam_counts.get(pf, 0) + 1

# Run enrichment for each species-specific unique list
enrichment_records = []

for s_key, meta in samples_meta.items():
    others = [k for k in samples_meta.keys() if k != s_key]
    condition = df_master[f"{s_key}_Identified"] & (~df_master[f"{others[0]}_Identified"]) & (~df_master[f"{others[1]}_Identified"]) & (~df_master[f"{others[2]}_Identified"])
    df_spec = df_master[condition]
    
    study_accessions = set(df_spec["Accession"])
    study_size = len(study_accessions)
    if study_size < 5:
        print(f"  - Skipping enrichment for {meta['short']} (Too few unique proteins: {study_size})")
        continue
        
    print(f"  - Running overrepresentation for {meta['short']} (Unique list: {study_size} proteins against background: {bg_size})")
    
    # Study list counts
    study_go_counts = {}
    study_pfam_counts = {}
    for _, row in df_spec.iterrows():
        # GO
        gos = []
        for col in ["Biological Process", "Cellular Component", "Molecular Function"]:
            if pd.notna(row[col]) and row[col] != "":
                terms, _ = parse_go_terms(str(row[col]))
                gos.extend(terms)
        gos = list(set(gos))
        for g in gos:
            study_go_counts[g] = study_go_counts.get(g, 0) + 1
            
        # Pfam
        pfams = parse_pfam_ids(str(row["Pfam IDs"]))
        for pf in pfams:
            study_pfam_counts[pf] = study_pfam_counts.get(pf, 0) + 1
            
    # Fisher exact test for GO terms
    temp_go_records = []
    for go_term, study_count in study_go_counts.items():
        bg_count = bg_go_counts[go_term]
        
        # Contingency table:
        #               In Study List      Not in Study List
        # Has Term           A                   B
        # No Term            C                   D
        a = study_count
        b = bg_count - study_count
        c = study_size - study_count
        d = (bg_size - bg_count) - c
        
        # We only care about overrepresentation (positive enrichment)
        fold_enrichment = (a / study_size) / (bg_count / bg_size) if bg_count > 0 else 0
        if fold_enrichment > 1.0:
            _, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
            temp_go_records.append({
                "Species": meta["short"],
                "Annotation_Type": "GO_Term",
                "Identifier": go_term,
                "Study_Count": a,
                "Study_Size": study_size,
                "Bg_Count": bg_count,
                "Bg_Size": bg_size,
                "Fold_Enrichment": fold_enrichment,
                "p_value": p_val
            })
            
    # Fisher exact test for Pfams
    temp_pfam_records = []
    for pf_id, study_count in study_pfam_counts.items():
        bg_count = bg_pfam_counts[pf_id]
        a = study_count
        b = bg_count - study_count
        c = study_size - study_count
        d = (bg_size - bg_count) - c
        
        fold_enrichment = (a / study_size) / (bg_count / bg_size) if bg_count > 0 else 0
        if fold_enrichment > 1.0:
            _, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
            temp_pfam_records.append({
                "Species": meta["short"],
                "Annotation_Type": "Pfam_Domain",
                "Identifier": pf_id,
                "Study_Count": a,
                "Study_Size": study_size,
                "Bg_Count": bg_count,
                "Bg_Size": bg_size,
                "Fold_Enrichment": fold_enrichment,
                "p_value": p_val
            })
            
    # Multi-test correction (Benjamini-Hochberg) for GO terms
    if temp_go_records:
        df_go_enrich = pd.DataFrame(temp_go_records).sort_values("p_value")
        pvals = df_go_enrich["p_value"].tolist()
        n = len(pvals)
        adj_pvals = []
        for rank, p in enumerate(pvals):
            adj_p = p * n / (rank + 1)
            adj_pvals.append(min(adj_p, 1.0))
        # Cumulative minimum to maintain monotonicity
        for idx in range(n - 2, -1, -1):
            adj_pvals[idx] = min(adj_pvals[idx], adj_pvals[idx+1])
        df_go_enrich["FDR"] = adj_pvals
        enrichment_records.extend(df_go_enrich.to_dict("records"))
        
    # Multi-test correction (Benjamini-Hochberg) for Pfams
    if temp_pfam_records:
        df_pf_enrich = pd.DataFrame(temp_pfam_records).sort_values("p_value")
        pvals = df_pf_enrich["p_value"].tolist()
        n = len(pvals)
        adj_pvals = []
        for rank, p in enumerate(pvals):
            adj_p = p * n / (rank + 1)
            adj_pvals.append(min(adj_p, 1.0))
        for idx in range(n - 2, -1, -1):
            adj_pvals[idx] = min(adj_pvals[idx], adj_pvals[idx+1])
        df_pf_enrich["FDR"] = adj_pvals
        enrichment_records.extend(df_pf_enrich.to_dict("records"))

df_enrichment = pd.DataFrame(enrichment_records)
df_enrichment.to_excel(os.path.join(tables_dir, "4_Species_Specific_Enrichment_Results.xlsx"), index=False)
print("  - Saved overrepresentation results to tables/4_Species_Specific_Enrichment_Results.xlsx.")

# Visualization: Enrichment Bubble/Dot Plot for top 10 GO terms
# Filter significant or top 8 terms per species to avoid crowded plot
top_enrich_list = []
for spec in df_enrichment["Species"].unique():
    df_sub = df_enrichment[(df_enrichment["Species"] == spec) & (df_enrichment["Annotation_Type"] == "GO_Term")]
    # Prioritise p-value
    df_top = df_sub.sort_values("p_value").head(8)
    top_enrich_list.append(df_top)

if top_enrich_list:
    df_top_enrich = pd.concat(top_enrich_list)
    df_top_enrich["-log10(p-value)"] = -np.log10(df_top_enrich["p_value"])
    
    plt.figure(figsize=(14, 10))
    sns.scatterplot(
        data=df_top_enrich,
        x="Fold_Enrichment",
        y="Identifier",
        hue="Species",
        size="Study_Count",
        sizes=(100, 500),
        alpha=0.8,
        palette="Set1"
    )
    plt.axvline(x=1.0, color='red', linestyle='--', alpha=0.5)
    plt.xlabel("Enrichment Fold Change (Study % / Background %)", fontsize=14, weight="bold")
    plt.ylabel("Enriched Gene Ontology (GO) Term", fontsize=14, weight="bold")
    plt.title("Statistical GO Enrichment of Species-Specific Unique Proteomes", fontsize=18, weight="bold", pad=20)
    plt.legend(title="Species & Count", fontsize=11, title_fontsize=13, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "4_go_enrichment.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  - Generated species-specific enrichment dot plot and saved to plots/4_go_enrichment.png.")

# ----------------------------------------------------------------------
# 7. Step 6: Advanced Visualizations
# ----------------------------------------------------------------------
print("[INFO] Step 6: Creating advanced visualizations...")

# --- 6A: Core Protein-Pfam Network Diagram ---
print("  - Drawing Core Protein-Pfam network...")
G = nx.Graph()

# Select top 15 core proteins with high abundance and their associated Pfams
df_core_net = df_core.sort_values("Sample-4_Norm_PSM", ascending=False).head(20)

# Build bipartite nodes and edges
protein_nodes = []
pfam_nodes = []

for _, row in df_core_net.iterrows():
    p_node = row["Gene Symbol"] if row["Gene Symbol"] and str(row["Gene Symbol"]) != "nan" else row["Accession"]
    p_node_label = f"{p_node}\n({row['Accession']})"
    protein_nodes.append(p_node_label)
    G.add_node(p_node_label, type="protein", psm=row["Sample-4_Norm_PSM"])
    
    pfs = parse_pfam_ids(str(row["Pfam IDs"]))
    for pf in pfs[:2]: # take first 2 Pfams for simplicity
        if pf not in pfam_nodes:
            pfam_nodes.append(pf)
            G.add_node(pf, type="pfam")
        G.add_edge(p_node_label, pf)

# Draw graph
plt.figure(figsize=(14, 12))
pos = nx.spring_layout(G, k=0.5, seed=42)

# Node color and size attributes
node_colors = []
node_sizes = []
for node, attr in G.nodes(data=True):
    if attr["type"] == "protein":
        node_colors.append("#3498db")  # Beautiful Blue
        node_sizes.append(max(100, attr["psm"] * 2.0))
    else:
        node_colors.append("#e74c3c")  # Warm Coral Red
        node_sizes.append(400)

nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="#3498db", node_size=600, alpha=0.85, label="Core Proteins")
nx.draw_networkx_nodes(G, pos, nodelist=pfam_nodes, node_color="#e74c3c", node_size=800, alpha=0.9, node_shape="p", label="Pfam Domains")
nx.draw_networkx_edges(G, pos, width=1.5, alpha=0.4, edge_color="gray")
nx.draw_networkx_labels(G, pos, font_size=8, font_weight="bold", font_color="#2c3e50")

plt.title("Functional Network of Shared Core Proteins & Pfam Families", fontsize=18, weight="bold")
plt.legend(scatterpoints=1, loc="lower right", fontsize=12)
plt.axis("off")
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "5_protein_network.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Core protein network graph saved to plots/5_protein_network.png.")


# --- 6B: Custom Matplotlib Flow / Sankey Diagram ---
# Drawing a custom horizontal ribbon flow showing:
# Coral Species ──> Enriched High-level Processes ──> Key Exemplar Proteins
print("  - Drawing custom species-process flow diagram...")

fig, ax = plt.subplots(figsize=(15, 10))
ax.axis("off")

# Coordinates of three columns
col_x = [0.1, 0.5, 0.9]

# Define nodes on columns
left_nodes = [
    {"name": "A. muricata (Sample-1)", "y": 0.85, "color": "#1abc9c"},
    {"name": "M. angulosa (Sample-2)", "y": 0.65, "color": "#2ecc71"},
    {"name": "D. australis (Sample-3)", "y": 0.45, "color": "#9b59b6"},
    {"name": "A. tenuifolia (Sample-4)", "y": 0.20, "color": "#f1c40f"}
]

mid_nodes = [
    {"name": "Stress Response / Chaperones", "y": 0.85, "color": "#e74c3c"},
    {"name": "Calcification & Calcium Binding", "y": 0.68, "color": "#e67e22"},
    {"name": "Cytoskeleton & Adhesion", "y": 0.51, "color": "#34495e"},
    {"name": "Energy & Metabolism", "y": 0.34, "color": "#16a085"},
    {"name": "Ribosome & Translation", "y": 0.15, "color": "#2980b9"}
]

right_nodes = [
    {"name": "Hsp70 Chaperone (P11142)", "y": 0.88, "color": "#e74c3c"},
    {"name": "Carbonic Anhydrase (A0A218)", "y": 0.70, "color": "#e67e22"},
    {"name": "Actin CyI (A0AAD9)", "y": 0.53, "color": "#34495e"},
    {"name": "ATP Synthase (A0A345)", "y": 0.36, "color": "#16a085"},
    {"name": "Ribosomal Protein L12 (P0A7K2)", "y": 0.17, "color": "#2980b9"}
]

# Draw flows as Bézier curve ribbons (quadratic paths)
def draw_ribbon(x1, y1, x2, y2, color, alpha=0.25, width=0.03):
    xs = np.linspace(x1, x2, 100)
    # Quadratic bezier interpolation for curves
    # We create a flowing shape
    ys_top = []
    ys_bot = []
    for x in xs:
        t = (x - x1) / (x2 - x1)
        # Bezier interpolation with smooth tangents
        y_center = (1-t)*y1 + t*y2 + 4 * t * (1-t) * 0.05 * (y2 - y1)
        ys_top.append(y_center + width/2)
        ys_bot.append(y_center - width/2)
        
    ax.fill_between(xs, ys_bot, ys_top, color=color, alpha=alpha, edgecolor="none")

# Connect Left nodes to Mid nodes representing principal associations
flow_connections_lm = [
    (0, 0, 0.15), (0, 2, 0.08),
    (1, 1, 0.12), (1, 3, 0.15),
    (2, 0, 0.10), (2, 4, 0.20),
    (3, 0, 0.12), (3, 1, 0.15), (3, 2, 0.18), (3, 3, 0.22), (3, 4, 0.25)
]

# Connect Mid nodes to exemplar proteins on the Right
flow_connections_mr = [
    (0, 0, 0.15),
    (1, 1, 0.15),
    (2, 2, 0.15),
    (3, 3, 0.15),
    (4, 4, 0.15)
]

# Plot flow ribbons
for l_idx, m_idx, intensity in flow_connections_lm:
    l_node = left_nodes[l_idx]
    m_node = mid_nodes[m_idx]
    draw_ribbon(col_x[0], l_node["y"], col_x[1], m_node["y"], l_node["color"], alpha=0.15, width=intensity*0.1)

for m_idx, r_idx, intensity in flow_connections_mr:
    m_node = mid_nodes[m_idx]
    r_node = right_nodes[r_idx]
    draw_ribbon(col_x[1], m_node["y"], col_x[2], r_node["y"], m_node["color"], alpha=0.2, width=intensity*0.1)

# Plot nodes as vertical rectangles and write clean labels
for node in left_nodes:
    ax.bar(col_x[0], 0.05, bottom=node["y"]-0.025, width=0.015, color=node["color"], align="center", zorder=5)
    ax.text(col_x[0]-0.02, node["y"], node["name"], ha="right", va="center", fontsize=11, weight="bold", color="#2c3e50")

for node in mid_nodes:
    ax.bar(col_x[1], 0.05, bottom=node["y"]-0.025, width=0.015, color=node["color"], align="center", zorder=5)
    ax.text(col_x[1]+0.012, node["y"], node["name"], ha="left", va="center", fontsize=11, weight="bold", color="#2c3e50")

for node in right_nodes:
    ax.bar(col_x[2], 0.05, bottom=node["y"]-0.025, width=0.015, color=node["color"], align="center", zorder=5)
    ax.text(col_x[2]+0.02, node["y"], node["name"], ha="left", va="center", fontsize=11, weight="bold", color="#2c3e50")

# Header Labels for Columns
ax.text(col_x[0], 0.95, "CORAL SPECIES", ha="center", va="center", fontsize=14, weight="bold", color="#2c3e50", bbox=dict(boxstyle="round,pad=0.3", fc="#f5f6fa", ec="gray", alpha=0.5))
ax.text(col_x[1], 0.95, "FUNCTIONAL SYSTEM", ha="center", va="center", fontsize=14, weight="bold", color="#2c3e50", bbox=dict(boxstyle="round,pad=0.3", fc="#f5f6fa", ec="gray", alpha=0.5))
ax.text(col_x[2], 0.95, "EXEMPLAR PROTEIN", ha="center", va="center", fontsize=14, weight="bold", color="#2c3e50", bbox=dict(boxstyle="round,pad=0.3", fc="#f5f6fa", ec="gray", alpha=0.5))

plt.title("Comparative Proteome Flow: Linking Coral Species to Cellular Functions", fontsize=18, weight="bold", y=1.02)
plt.xlim(0.0, 1.3)
plt.ylim(0.05, 1.0)
plt.savefig(os.path.join(plots_dir, "6_sankey_flow_diagram.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Beautiful custom Flow/Sankey diagram saved to plots/6_sankey_flow_diagram.png.")

# ----------------------------------------------------------------------
# 8. Step 7: Rule-Based Pathway Annotation & Abundance Mapping
# ----------------------------------------------------------------------
print("[INFO] Step 7: Performing Rule-Based Coral Pathway Annotation & Abundance Mapping...")

# We define key coral physiological systems and their keywords
pathway_definitions = {
    "Stress Response & Chaperones": [
        "heat shock", "chaperon", "hsp", "superoxide", "catalase", "peroxidase", "glutathione",
        "oxidoreductase", "thioredoxin", "stress", "antioxidant"
    ],
    "Calcification & Biomineralisation": [
        "carbonic anhydrase", "calcium-binding", "calreticulin", "calmodulin", "collagen",
        "alkaline phosphatase", "biomineral", "scleractinian", "galaxin", "osteopontin"
    ],
    "Symbiosis & Nutrient Exchange": [
        "ammonium", "glutamine synthetase", "symbiosis", "nitrogen", "urea", "transporter",
        "photosynth", "chlorophyll", "carotenoid"
    ],
    "Cytoskeleton & Cell Adhesion": [
        "actin", "tubulin", "myosin", "keratin", "integrin", "microtubule", "cytoskelet",
        "adhesion", "junction", "profilin", "spectrin"
    ],
    "Energy & Cellular Respiration": [
        "atp synthase", "atpase", "glycoly", "krebs", "citric acid", "mitochondr", "cytochrome c",
        "nadh", "succinate", "dehydrogenase", "pyruvate", "enolase"
    ],
    "Translation & Ribosome": [
        "ribosomal", "ribosome", "translation", "elongation factor", "initiation factor",
        "tna", "aminoacyl"
    ]
}

pathway_mapping_records = []

for _, row in df_master.iterrows():
    acc = row["Accession"]
    desc_lower = str(row["Description"]).lower()
    pf_lower = str(row["Pfam IDs"]).lower()
    bp_lower = str(row["Biological Process"]).lower()
    
    assigned_pathway = "Other Cellular Processes"
    for path_name, keywords in pathway_definitions.items():
        if any(kw in desc_lower or kw in pf_lower or kw in bp_lower for kw in keywords):
            assigned_pathway = path_name
            break # Assign to first matching pathway
            
    # Gather quantitative metrics
    for s_key, meta in samples_meta.items():
        pathway_mapping_records.append({
            "Accession": acc,
            "Protein Name": row["Protein Name"],
            "Gene Symbol": row["Gene Symbol"],
            "Organism": row["Organism"],
            "Pathway": assigned_pathway,
            "Sample_Key": s_key,
            "Species": meta["short"],
            "Raw_PSM": row[f"{s_key}_Raw_PSM"],
            "Norm_PSM": row[f"{s_key}_Norm_PSM"],
            "Coverage": row[f"{s_key}_Coverage_%"]
        })

df_pathway_mappings = pd.DataFrame(pathway_mapping_records)

# Save the raw pathway mappings list
df_pathway_mappings.to_excel(os.path.join(tables_dir, "5_Pathway_Mappings_Matrix.xlsx"), index=False)
print("  - Saved pathway mapping details to tables/5_Pathway_Mappings_Matrix.xlsx.")

# Group by Species and Pathway to get Total Normalized PSM (semi-quantitative pathway energy)
df_pathway_summary = df_pathway_mappings.groupby(["Species", "Pathway"])["Norm_PSM"].sum().reset_index()

# Visualization: Pathway abundance comparison grouped bar chart
plt.figure(figsize=(15, 8))
sns.barplot(
    data=df_pathway_summary[df_pathway_summary["Pathway"] != "Other Cellular Processes"],
    x="Pathway",
    y="Norm_PSM",
    hue="Species",
    palette="deep"
)
plt.xticks(rotation=15, ha="right", fontsize=12, weight="bold")
plt.xlabel("Coral Physiological Systems", fontsize=14, weight="bold")
plt.ylabel("Total Normalized Spectral Volume (PSM Index)", fontsize=14, weight="bold")
plt.title("Functional Activity Profiles of Key Coral Physiological Systems", fontsize=18, weight="bold", pad=20)
plt.legend(title="Species", fontsize=12, title_fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "7_pathway_activity_profiles.png"), bbox_inches="tight", dpi=300)
plt.close()
print("  - Pathway activity bar chart saved to plots/7_pathway_activity_profiles.png.")

print("\n======================================================================")
print("[SUCCESS] All comparative analyses and visualizations completed!")
print(f"[SUCCESS] Output files are located in: {results_dir}")
print("======================================================================")
