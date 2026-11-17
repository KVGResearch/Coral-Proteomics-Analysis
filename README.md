# Coral Proteomics Analysis

A reproducible Python workflow for comparative analysis of coral tissue proteomics datasets from the **Gulf of Mannar, Tamil Nadu, India**. The pipeline processes protein-identification tables, compares shared and species-specific proteins, retrieves functional annotations, performs exploratory enrichment analyses, and generates figures with accompanying CSV data.

## Study samples

| Sample | Species | Expected input file |
|---|---|---|
| S1 | *Acropora muricata* | `Veena_Sample-1_protein_list.xlsx` |
| S2 | *Mussa angulosa* | `Veena_Sample-2_protein_list.xlsx` |
| S3 | *Dendronephthya australis* | `Veena_Sample-3_protein_list.xlsx` |
| S4 | *Agaricia tenuifolia* | `Veena_Sample-4_protein_list.xlsx` |

The inputs are Excel protein lists exported from Proteome Discoverer. The script retains proteins marked `High` in `Protein FDR Confidence: Combined` and excludes entries marked as contaminants when those columns exist.

## Analysis features

- Protein counts, sequence coverage, peptide counts, molecular weight, and calculated isoelectric point
- Shared, unique, and four-species core protein comparisons, including UpSet plots and Jaccard similarity
- Gene Ontology (GO) Biological Process, Cellular Component, and Molecular Function annotations
- Reactome, KEGG, and optionally WikiPathways membership and enrichment analyses
- Comparative visualizations: heatmaps, Sankey-style plots, chord diagrams, and networks
- Export of numerical CSV data underlying the generated figures

**Note:** Network figures use shared annotation/pathway information to infer co-associations; they should not be interpreted as experimentally validated protein–protein interactions.

## Requirements

- Python 3
- Internet access for first-time annotation retrieval
- Four Excel input files with the expected Proteome Discoverer columns

Install Python dependencies:

```bash
pip install numpy pandas openpyxl matplotlib seaborn networkx scipy statsmodels
```

## Usage

1. Place the analysis script and four Excel files in the same directory.
2. Update the `SAMPLE_FILES` mapping if your Excel filenames differ.
3. Run the pipeline:

```bash
python coral_proteomics_full_pipeline.py
```

The script expects fields including `Accession`, `Description`, `Coverage [%]`, `# Peptides`, `MW [kDa]`, and `calc. pI`. Check the column-name configuration in the script if your exports differ.

## Outputs

The current script writes results to:

```text
coral_figures/
├── Fig01_*.png ... Fig31_*.png
└── csv_data/
    ├── Fig01_*.csv ...
    ├── Top_GO_Terms_Summary.csv
    └── Annotation_Fetch_Coverage.csv
annotation_cache.json
```

It is configured to produce **31 figures** at 180 dpi, together with supporting CSV files. Figures cover protein identification, comparative proteomics, GO classification, pathway enrichment, core and unique proteins, and annotation-based network visualizations.

## Annotation and statistical methods

The pipeline uses UniProt accessions to retrieve GO terms and links to KEGG and Reactome. WikiPathways querying is **disabled by default** (`FETCH_WIKIPATHWAYS = False`) because individual requests can be slow or inaccessible. Set it to `True` if required and network access permits.

Annotations are cached in `annotation_cache.json` to reduce repeated requests. Enrichment uses one-sided Fisher's exact tests and Benjamini–Hochberg multiple-testing correction, with the pooled four-species protein set as background.

**Interpretation caution:** The current code uses exploratory adjusted-p-value cutoffs (including 0.5 and 0.9 in some analyses), which are much less stringent than conventional FDR thresholds. Review these settings and the background universe before interpreting results as statistically significant. Missing database annotations do not imply absence of biological function.

## Authors

**Santhosh M E**  
Doctoral Researcher  
Vellore Institute of Technology (VIT), Vellore  
Email: santhosh.me@vit.ac.in

**Dr. George Priya Doss C**  
Professor  
Vellore Institute of Technology (VIT), Vellore  
Email: georgepriyadoss@vit.ac.in

## Data sharing

Check project permissions before publishing sample-level data or raw mass spectrometry files. Large files can be stored outside GitHub, with access instructions added here where appropriate.

## License

The analysis code is intended for release under the **MIT License**. See the repository `LICENSE` file once added. This license applies to the code, not automatically to experimental datasets, third-party databases, or manuscript figures. Confirm the copyright holder(s) and institutional permissions before public release.

## Citation

If you use this workflow, please cite the associated research article once published. The article citation and DOI will be added here after publication. For reproducibility, cite the specific repository release or archived DOI used in your analysis.
