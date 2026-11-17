#!/usr/bin/env python3
"""
=============================================================================
  CORAL TISSUE PROTEOMICS – COMPLETE ANALYSIS PIPELINE (Figures 1–31)
  Author 
  1. Santhosh M E
     Doctoral Researcher
     (santhosh.me@vit.ac.in)
  2. Dr George Priya Doss C
     Professor
     (georgepriyadoss@vit.ac.in)
=============================================================================

  USAGE
  -----
  1. Place this script in the SAME FOLDER as your four xlsx files.
  2. Edit SAMPLE_FILES below to match your exact filenames.
  3. Install dependencies:
       pip install pandas openpyxl matplotlib seaborn networkx scipy statsmodels
  4. Run:
       python coral_proteomics_full_pipeline.py
  5. INTERNET REQUIRED on first run. Every protein's real GO terms and
     Reactome/KEGG/WikiPathways memberships are fetched live from UniProt,
     KEGG and Reactome/WikiPathways web services using the "Accession"
     column in your Excel files (expects UniProt-style accessions, e.g.
     P12345 / A0A2B4XYZ1). Results are cached to annotation_cache.json in
     this folder, so re-runs after the first one are fast and work offline.

  OUTPUT
  ------
  Folder "coral_figures/" with 31 publication-quality PNG figures (180 dpi).
  Folder "coral_figures/csv_data/" with one CSV per figure containing the
  exact numbers plotted, so every figure can be cross-checked against data.
  File "annotation_cache.json" — cached real annotations per accession, so
  you are not re-querying the same protein IDs on every run.

  FIGURES GENERATED
  -----------------
  Fig 01  Per-sample protein summary (MW vs Coverage scatter, 2x2)
  Fig 02  Comparative metrics (count bar, coverage violin, peptide box)
  Fig 03  UpSet plot (all intersection combinations)
  Fig 04  Jaccard similarity heatmap
  Fig 05  GO-BP stacked bar per species
  Fig 06  GO-CC and GO-MF bubble charts
  Fig 07  GO comparative heatmap (all 3 ontologies)
  Fig 08  Sankey / alluvial: species → biological process
  Fig 09  Chord diagram: protein sharing between species
  Fig 10  Bipartite functional network: species × GO-BP terms
  Fig 11  pI vs MW scatter (all species overlaid)
  Fig 12  MW and pI overlay histograms
  Fig 13  Pathway keyword analysis (grouped bar)
  Fig 14  Cross-species coverage heatmap (top 40 shared proteins)
  Fig 15  Radar chart: normalised proteomic profile
  Fig 16  GO-BP pie charts per species (2x2)
  Fig 17  Species-unique proteins top-15 (horizontal bars)
  Fig 18  Comprehensive summary statistics table
  Fig 19  Full analysis workflow diagram
  Fig 20  Reactome enrichment per species (Fisher's exact, BH-FDR)
  Fig 21  WikiPathways enrichment per species
  Fig 22  KEGG enrichment per species
  Fig 23  Comparative Reactome heatmap (all 4 species)
  Fig 24  Comparative KEGG heatmap
  Fig 25  Comparative WikiPathways heatmap
  Fig 26  Species-specific GO enrichment of unique proteins (4x3)
  Fig 27  Species-specific pathway enrichment of unique proteins (4x3)
  Fig 28  Core proteome GO + pathway enrichment summary
  Fig 29  PPI co-association network – core proteome
  Fig 30  Chord diagram: species × enriched pathways
  Fig 31  STRING-style PPI network – top 40 biologically important proteins

  NOTE ON ENRICHMENT
  ------------------
  GO terms (BP/CC/MF) and Reactome/KEGG/WikiPathways pathway memberships are
  fetched live per-protein from UniProt (rest.uniprot.org), KEGG
  (rest.kegg.jp), Reactome (reactome.org/ContentService) and WikiPathways
  (webservice.wikipathways.org), keyed on the "Accession" column of your
  input files. Nothing is keyword-guessed from free-text descriptions.
  Fisher's exact test + Benjamini-Hochberg FDR is then run per species
  (study set) against the pooled 4-species background, over whatever real
  pathway/GO terms were actually returned by these databases. Coral proteins
  with no curated hit in a given database simply contribute no membership
  for it — that's a real (if sparse) result, not a workaround.

=============================================================================
"""

import os, sys, math, json, time, re, warnings
import urllib.request, urllib.parse, urllib.error
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import networkx as nx
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
#  ★  FILE CONFIGURATION  – edit filenames here  ★
# ─────────────────────────────────────────────────────────────────────────────
SAMPLE_FILES = {
    "S1": {"species": "Acropora muricata",        "color": "#1B4F72",
           "file": "Veena_Sample-1_protein_list.xlsx"},
    "S2": {"species": "Mussa angulosa",           "color": "#1E8449",
           "file": "Veena_Sample-2_protein_list.xlsx"},
    "S3": {"species": "Dendronephthya australis", "color": "#784212",
           "file": "Veena_Sample-3_protein_list.xlsx"},
    "S4": {"species": "Agaricia tenuifolia",      "color": "#7D3C98",
           "file": "Veena_Sample-4_protein_list.xlsx"},
}

OUT_DIR = "coral_figures_17-07_V2"
CSV_DIR = os.path.join(OUT_DIR, "csv_data")

# WikiPathways can only be queried ONE protein at a time (no reliable batch
# lookup), so for a few thousand proteins it can take a very long time and is
# the most likely to get stuck behind an institutional proxy. It's OFF by
# default — GO (UniProt) + Reactome + KEGG enrichment all still run fully.
# Set to True only if you have a fast, unrestricted connection and are happy
# to wait (roughly 3 proteins/second).
FETCH_WIKIPATHWAYS = False
TOP_N_TERMS = 15   # how many top BP/CC/MF terms to report in the summary
TOP_N_GO_PLOT = 15 # how many top GO terms each chart (Fig 5/6/7/8/16) shows before bucketing the rest as "Other"

# Proteome Discoverer column names (change only if your export differs)
COL_ACCESSION   = "Accession"
COL_DESC        = "Description"
COL_COVERAGE    = "Coverage [%]"
COL_PEPTIDES    = "# Peptides"
COL_PSM         = "# PSMs"
COL_UNIQ_PEP   = "# Unique Peptides"
COL_AAS         = "# AAs"
COL_MW          = "MW [kDa]"
COL_PI          = "calc. pI"
COL_PEP_SCORE   = "Sum PEP Score"
COL_BP          = "Biological Process"
COL_CC          = "Cellular Component"
COL_MF          = "Molecular Function"
COL_FDR         = "Protein FDR Confidence: Combined"
COL_CONTAMINANT = "Contaminant"

# ─────────────────────────────────────────────────────────────────────────────
#  LIVE ANNOTATION FETCHING  (UniProt · KEGG · Reactome · WikiPathways)
# ─────────────────────────────────────────────────────────────────────────────
CACHE_FILE       = "annotation_cache.json"
UNIPROT_STREAM   = "https://rest.uniprot.org/uniprotkb/stream"
KEGG_BASE        = "https://rest.kegg.jp"
REACTOME_BASE    = "https://reactome.org/ContentService"
WIKIPATHWAYS_BASE= "https://webservice.wikipathways.org"
HTTP_HEADERS     = {"User-Agent": "CoralProteomicsPipeline/1.0 (research script)"}
UNIPROT_CHUNK    = 90     # accessions per UniProt query (stays well under URL length limits)
KEGG_CHUNK       = 10     # gene ids per KEGG link/pathway request (KEGG fair-use guidance)
REACTOME_CHUNK   = 40     # ids per Reactome batch POST
REQUEST_SLEEP    = 0.34   # ~3 requests/sec, polite to all four services

def _http_get(url, timeout=20, retries=2, diagnose=False):
    """Returns response text, or None on failure. Permanent errors (4xx) are NOT retried —
    only transient network errors/timeouts are. Pass diagnose=True to print status + a body
    snippet even on failure, useful for figuring out *why* a request is failing."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HTTP_HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
                if diagnose:
                    print(f"      [DEBUG] GET {url[:110]}\n              -> HTTP {resp.status}, "
                          f"{len(text)} bytes, starts with: {text[:150]!r}")
                return text
        except urllib.error.HTTPError as e:
            body = ""
            try: body = e.read().decode("utf-8", errors="replace")[:150]
            except Exception: pass
            print(f"      [WARN] GET {url[:90]}... -> HTTP {e.code} {e.reason} "
                  f"(not retrying, this is a permanent error): {body!r}")
            return None  # 4xx/5xx from the server itself — retrying won't help
        except Exception as e:
            if attempt == retries - 1:
                print(f"      [WARN] GET failed ({url[:90]}...): {e}")
                return None
            time.sleep(1.5 * (attempt + 1))
    return None

def _http_post(url, body_text, timeout=20, retries=2):
    data = body_text.encode("utf-8")
    headers = dict(HTTP_HEADERS); headers.update({"Content-Type": "text/plain", "Accept": "application/json"})
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            body = ""
            try: body = e.read().decode("utf-8", errors="replace")[:150]
            except Exception: pass
            print(f"      [WARN] POST {url[:90]}... -> HTTP {e.code} {e.reason} "
                  f"(not retrying): {body!r}")
            return None
        except Exception as e:
            if attempt == retries - 1:
                print(f"      [WARN] POST failed ({url[:90]}...): {e}")
                return None
            time.sleep(1.5 * (attempt + 1))
    return None

def load_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, "r") as f: return json.load(f)
        except Exception:
            return {}
    return {}

def save_cache(cache):
    try:
        with open(CACHE_FILE, "w") as f: json.dump(cache, f)
    except Exception as e:
        print(f"      [WARN] could not write {CACHE_FILE}: {e}")

def _chunks(seq, n):
    for i in range(0, len(seq), n): yield seq[i:i+n]

def _clean_go(field_val):
    """UniProt tsv go_p/go_c/go_f cells look like 'term one [GO:0006914]; term two [GO:0007005]'."""
    out = []
    if not field_val: return out
    for seg in field_val.split(";"):
        seg = seg.strip()
        if not seg: continue
        m = re.match(r"^(.*?)\s*\[GO:\d+\]$", seg)
        out.append(m.group(1).strip() if m else seg)
    return out

def fetch_uniprot_batch(accessions, cache):
    """Populate cache[acc] = {go_bp, go_cc, go_mf, kegg_genes, reactome_ids} for accessions not yet cached."""
    todo = [a for a in accessions if a not in cache or "go_bp" not in cache.get(a, {})]
    if not todo:
        return
    print(f"      Fetching UniProt annotations for {len(todo)} accessions ({UNIPROT_CHUNK}/request)…")
    found = 0
    chunk_list = list(_chunks(todo, UNIPROT_CHUNK))
    for ci, chunk in enumerate(chunk_list):
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = (f"{UNIPROT_STREAM}?query={urllib.parse.quote(query)}"
               f"&fields=accession,go_p,go_c,go_f,xref_kegg,xref_reactome&format=tsv")
        text = _http_get(url, diagnose=(ci == 0))  # print a raw snippet for the very first request
        if text:
            lines = text.strip().split("\n")
            for line in lines[1:]:  # skip header
                parts = line.split("\t")
                if len(parts) < 6: continue
                acc, go_p, go_c, go_f, xk, xr = parts[:6]
                cache[acc] = {
                    "go_bp": _clean_go(go_p), "go_cc": _clean_go(go_c), "go_mf": _clean_go(go_f),
                    "kegg_genes": [g for g in xk.split(";") if g.strip()],
                    "reactome_ids": [r for r in xr.split(";") if r.strip()],
                }
                found += 1
        if (ci + 1) % 10 == 0: print(f"        …{ci+1}/{len(chunk_list)} UniProt requests done")
        # accessions UniProt didn't return anything for (obsolete/invalid ID, or truly no data) → empty record
        for a in chunk:
            if a not in cache:
                cache[a] = {"go_bp": [], "go_cc": [], "go_mf": [], "kegg_genes": [], "reactome_ids": []}
        time.sleep(REQUEST_SLEEP)
    print(f"      ↳ {found}/{len(todo)} accessions returned UniProt data")

def fetch_kegg_pathway_names(cache):
    """Map each accession's kegg_genes -> real KEGG pathway (id, name) pairs."""
    gene_ids = sorted({g for rec in cache.values() for g in rec.get("kegg_genes", [])})
    if not gene_ids:
        return
    need_link = any("kegg_pathways" not in rec for rec in cache.values() if rec.get("kegg_genes"))
    if not need_link:
        return
    print(f"      Fetching KEGG pathway links for {len(gene_ids)} genes…")
    gene_to_paths = {}
    for chunk in _chunks(gene_ids, KEGG_CHUNK):
        url = f"{KEGG_BASE}/link/pathway/{'+'.join(chunk)}"
        text = _http_get(url)
        if text:
            for line in text.strip().split("\n"):
                if not line.strip() or "\t" not in line: continue
                gene, path_id = line.split("\t")
                gene_to_paths.setdefault(gene, []).append(path_id.replace("path:", ""))
        time.sleep(REQUEST_SLEEP)
    # Generic pathway-name dictionary (one call covers ~550 "map#####" pathways)
    name_map = {}
    listing = _http_get(f"{KEGG_BASE}/list/pathway")
    if listing:
        for line in listing.strip().split("\n"):
            if "\t" not in line: continue
            pid, pname_ = line.split("\t")
            name_map[pid.replace("path:", "")] = pname_.strip()
    def resolve_name(pid):
        if pid in name_map: return name_map[pid]
        m = re.search(r"(\d{5})$", pid)
        if m: return name_map.get("map" + m.group(1), pid)
        return pid
    for acc, rec in cache.items():
        if "kegg_pathways" in rec: continue
        pairs = []
        for g in rec.get("kegg_genes", []):
            for pid in gene_to_paths.get(g, []):
                pairs.append([pid, resolve_name(pid)])
        rec["kegg_pathways"] = pairs

def fetch_reactome_pathway_names(cache):
    """Resolve every distinct Reactome stable ID -> displayName via one batched POST call per chunk."""
    ids = sorted({r for rec in cache.values() for r in rec.get("reactome_ids", [])})
    if not ids: return
    need = any("reactome_pathways" not in rec for rec in cache.values() if rec.get("reactome_ids"))
    if not need: return
    print(f"      Fetching Reactome pathway names for {len(ids)} IDs…")
    id_to_name = {}
    for chunk in _chunks(ids, REACTOME_CHUNK):
        text = _http_post(f"{REACTOME_BASE}/data/query/ids", ",".join(chunk))
        if text:
            try:
                data = json.loads(text)
                if isinstance(data, dict): data = [data]
                for entry in data:
                    rid = entry.get("stId") or entry.get("dbId")
                    nm  = entry.get("displayName")
                    if rid and nm: id_to_name[str(rid)] = nm
            except Exception as e:
                print(f"      [WARN] could not parse Reactome response: {e}")
        time.sleep(REQUEST_SLEEP)
    for acc, rec in cache.items():
        if "reactome_pathways" in rec: continue
        rec["reactome_pathways"] = [[rid, id_to_name.get(rid, rid)] for rid in rec.get("reactome_ids", [])]

def fetch_wikipathways(accessions, cache):
    """Per-accession lookup (WikiPathways doesn't reliably map multi-id batch results back to a source id)."""
    todo = [a for a in accessions if "wikipathways" not in cache.get(a, {})]
    if not todo: return
    print(f"      Fetching WikiPathways memberships for {len(todo)} accessions (this is the slow one)…")
    consecutive_failures = 0
    for i, acc in enumerate(todo, 1):
        url = f"{WIKIPATHWAYS_BASE}/findPathwaysByXref?ids={acc}&codes=S&format=json"
        text = _http_get(url, timeout=12, retries=1)
        pairs = []
        if text:
            consecutive_failures = 0
            try:
                data = json.loads(text)
                for r in data.get("result", []):
                    if r.get("id") and r.get("name"): pairs.append([r["id"], r["name"]])
            except Exception:
                pass
        else:
            consecutive_failures += 1
        cache.setdefault(acc, {})["wikipathways"] = pairs
        if i % 50 == 0: print(f"        …{i}/{len(todo)}")
        if consecutive_failures >= 15:
            print(f"      [WARNING] WikiPathways failed {consecutive_failures} times in a row — "
                  f"the endpoint looks unreachable from this network. Stopping WikiPathways fetch "
                  f"early ({i}/{len(todo)} attempted); everything else continues normally.")
            for a in todo[i:]: cache.setdefault(a, {}).setdefault("wikipathways", [])
            break
        time.sleep(REQUEST_SLEEP)

def top_go_terms_summary(dfs, keys, shorts, top_n=TOP_N_TERMS):
    """Prints + saves the top N BP/CC/MF terms overall and per species. Called right after the
    fetch step so you have real, checkable output immediately without waiting for all 31 figures."""
    print(f"\n      TOP {top_n} GO TERMS (pooled across all species)")
    csv_rows = []
    for col, label in [(COL_BP,"Biological Process"),(COL_CC,"Cellular Component"),(COL_MF,"Molecular Function")]:
        pooled = Counter()
        for k in keys: pooled.update(parse_go(dfs[k][col]))
        top = pooled.most_common(top_n)
        print(f"      --- {label} ---")
        for term, cnt in top: print(f"        {cnt:>4}  {term}")
        for term, cnt in top: csv_rows.append({"ontology":label,"scope":"ALL_SPECIES","term":term,"protein_count":cnt})
        for si,k in enumerate(keys):
            per_sp = parse_go(dfs[k][col]).most_common(top_n)
            for term, cnt in per_sp:
                csv_rows.append({"ontology":label,"scope":shorts[si],"term":term,"protein_count":cnt})
    save_csv("Top_GO_Terms_Summary.csv", pd.DataFrame(csv_rows))
    print(f"\n      ↳ full per-species top-{top_n} tables saved to csv_data/Top_GO_Terms_Summary.csv\n")

def build_real_annotations(dfs):
    """Fetch (or load from cache) real GO + pathway annotations for every accession in the dataset,
    overwrite each df's GO columns with the real terms, and return per-database membership maps
    (accession -> set of "Name (ID)" strings) used by the enrichment engine."""
    all_accs = sorted(set().union(*[set(df[COL_ACCESSION]) for df in dfs.values()]))
    cache = load_cache()
    fetch_uniprot_batch(all_accs, cache)
    fetch_kegg_pathway_names(cache)
    fetch_reactome_pathway_names(cache)
    if FETCH_WIKIPATHWAYS:
        fetch_wikipathways(all_accs, cache)
    else:
        print("      Skipping WikiPathways (FETCH_WIKIPATHWAYS=False — see config at top of file).")
        for a in all_accs:
            cache.setdefault(a, {}).setdefault("wikipathways", [])
    save_cache(cache)

    n_with_go   = sum(1 for a in all_accs if cache.get(a, {}).get("go_bp") or cache.get(a, {}).get("go_cc") or cache.get(a, {}).get("go_mf"))
    n_with_kegg = sum(1 for a in all_accs if cache.get(a, {}).get("kegg_pathways"))
    n_with_reac = sum(1 for a in all_accs if cache.get(a, {}).get("reactome_pathways"))
    n_with_wp   = sum(1 for a in all_accs if cache.get(a, {}).get("wikipathways"))
    print(f"      Annotation coverage: GO {n_with_go}/{len(all_accs)} | "
          f"KEGG {n_with_kegg}/{len(all_accs)} | Reactome {n_with_reac}/{len(all_accs)} | "
          f"WikiPathways {n_with_wp}/{len(all_accs)}")
    if n_with_go == 0 and n_with_kegg == 0 and n_with_reac == 0:
        print("\n      [WARNING] No GO/KEGG/Reactome annotations were retrieved — likely your network "
              "blocks or proxies\n      rest.uniprot.org / rest.kegg.jp / reactome.org. The 31 figures will "
              "still be generated,\n      but pathway/GO-based ones (5,6,7,13,16,20-28,29,31) will show "
              "'no data'.\n      Test in a browser: "
              "https://rest.uniprot.org/uniprotkb/search?query=accession:P05067&fields=accession,go_p&format=tsv\n"
              "      If that page doesn't load, this machine's network is blocking it — try a different "
              "network,\n      or run this script on a machine with unrestricted internet.\n")

    for k, df in dfs.items():
        df[COL_BP] = df[COL_ACCESSION].map(lambda a: ";".join(cache.get(a, {}).get("go_bp", [])))
        df[COL_CC] = df[COL_ACCESSION].map(lambda a: ";".join(cache.get(a, {}).get("go_cc", [])))
        df[COL_MF] = df[COL_ACCESSION].map(lambda a: ";".join(cache.get(a, {}).get("go_mf", [])))

    def membership_map(field):
        m = {}
        for a in all_accs:
            pairs = cache.get(a, {}).get(field, [])
            m[a] = {f"{name} ({pid})" for pid, name in pairs}
        return m
    reactome_mem = membership_map("reactome_pathways")
    kegg_mem     = membership_map("kegg_pathways")
    wp_mem       = membership_map("wikipathways")
    coverage = {"accession": all_accs,
                "has_GO": [bool(cache.get(a,{}).get("go_bp") or cache.get(a,{}).get("go_cc") or cache.get(a,{}).get("go_mf")) for a in all_accs],
                "n_KEGG_pathways": [len(kegg_mem[a]) for a in all_accs],
                "n_Reactome_pathways": [len(reactome_mem[a]) for a in all_accs],
                "n_WikiPathways": [len(wp_mem[a]) for a in all_accs]}
    return dfs, reactome_mem, kegg_mem, wp_mem, coverage

# ─────────────────────────────────────────────────────────────────────────────
#  HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def sn(sp):
    p = sp.split(); return f"{p[0][0]}. {' '.join(p[1:])}"

def parse_go(series):
    c = Counter()
    for val in series.dropna():
        for t in str(val).split(";"):
            t = t.strip()
            if t and t.lower() != "nan": c[t] += 1
    return c

def top_terms_with_other(pooled_counts, n=TOP_N_GO_PLOT, other_label="Other terms"):
    """Given a Counter/dict of term->pooled count, keep only the top n terms; everything else
    is collapsed into one 'Other (k terms)' bucket so charts stay readable with real GO data
    (which can easily contain hundreds of distinct terms)."""
    ranked = sorted(pooled_counts.items(), key=lambda x: -x[1])
    keep = [t for t, _ in ranked[:n]]
    n_other = len(ranked) - len(keep)
    other_key = f"{other_label} ({n_other})" if n_other > 0 else None
    return keep, other_key

def short_desc(desc, n=40):
    d = str(desc).split("OS=")[0].strip()
    return d[:n] + "…" if len(d) > n else d

def pname(desc):
    return str(desc).split(" OS=")[0].strip()

def pw_short(name, n=50):
    return name.split(" (")[0][:n]

def savefig(name, fig):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {name}")

def save_csv(name, df):
    """Write the exact data behind a figure to coral_figures/csv_data/<name>
    so every plot can be cross-checked against numbers."""
    path = os.path.join(CSV_DIR, name)
    try:
        if isinstance(df, pd.DataFrame):
            df.to_csv(path, index=False)
        else:
            pd.DataFrame(df).to_csv(path, index=False)
        print(f"      ↳ csv: {name}")
    except Exception as e:
        print(f"      [WARN] could not write {name}: {e}")

def get_go_set(row):
    g = set()
    for col in [COL_BP, COL_CC, COL_MF]:
        for t in str(row.get(col, "")).split(";"):
            t = t.strip()
            if t and t.lower() != "nan": g.add(t)
    return g

# ─────────────────────────────────────────────────────────────────────────────
#  LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────

def load_samples():
    dfs = {}
    for key, m in SAMPLE_FILES.items():
        if not os.path.exists(m["file"]):
            sys.exit(f"\n[ERROR] File not found: {m['file']}\n"
                     f"        Edit SAMPLE_FILES at the top of this script.\n")
        df = pd.read_excel(m["file"])
        if COL_FDR in df.columns:
            df = df[df[COL_FDR] == "High"]
        if COL_CONTAMINANT in df.columns:
            df = df[df[COL_CONTAMINANT] == False]
        dfs[key] = df.reset_index(drop=True)
        print(f"  Loaded {key} ({m['species']}): {len(dfs[key])} proteins")
    return dfs

# ─────────────────────────────────────────────────────────────────────────────
#  ENRICHMENT ENGINE  (Fisher's exact + BH-FDR over REAL fetched memberships)
# ─────────────────────────────────────────────────────────────────────────────

def pathway_fisher_enrichment(df_study, df_bg, membership, fdr_thr=0.5):
    """membership: accession -> set of real pathway/term names (already fetched from a live DB)."""
    n_s = len(df_study); n_b = len(df_bg)
    study_accs = list(df_study[COL_ACCESSION]); bg_accs = list(df_bg[COL_ACCESSION])
    universe = set()
    for a in bg_accs: universe |= membership.get(a, set())
    rows = []
    for term in universe:
        ks = sum(1 for a in study_accs if term in membership.get(a, set()))
        kb = sum(1 for a in bg_accs if term in membership.get(a, set()))
        if ks == 0: continue
        a, b, c, d = ks, n_s-ks, kb-ks, n_b-n_s-(kb-ks)
        if d < 0: continue
        try: _, pv = fisher_exact([[a,b],[c,d]], alternative="greater")
        except: pv = 1.0
        rows.append({"pathway": term, "study_hits": a, "bg_hits": kb,
                     "study_n": n_s, "bg_n": n_b, "pval": pv,
                     "ratio": a/n_s if n_s > 0 else 0})
    if not rows: return pd.DataFrame()
    res = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(res["pval"], method="fdr_bh") if len(res)>1 else ([None],[res["pval"].values[0]],[None],[None])
    res["padj"] = padj
    return res[res["padj"] < fdr_thr].sort_values("padj")

def go_fisher_enrichment(df_study, df_full, go_col, fdr_thr=0.9):
    n_s = len(df_study); n_f = len(df_full)
    uc = parse_go(df_study[go_col]); fc = parse_go(df_full[go_col])
    rows = []
    for term, ku in uc.items():
        kf = fc.get(term, 0)
        a, b, c, d = ku, n_s-ku, kf-ku, n_f-n_s-(kf-ku)
        if a == 0 or d < 0: continue
        try: _, pv = fisher_exact([[a,b],[c,d]], alternative="greater")
        except: pv = 1.0
        rows.append({"term": term, "uniq_hits": ku, "full_hits": kf, "pval": pv})
    if not rows: return pd.DataFrame()
    res = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(res["pval"], method="fdr_bh") if len(res)>1 else ([None],[res["pval"].values[0]],[None],[None])
    res["padj"] = padj
    return res[res["padj"] < fdr_thr].sort_values("padj")

def run_all_enrichment(dfs, keys, reactome_mem, kegg_mem, wp_mem):
    all_dbs = {"REACTOME": reactome_mem, "KEGG": kegg_mem, "WIKIPATHWAYS": wp_mem}
    bg = pd.concat(list(dfs.values()), ignore_index=True).drop_duplicates(subset=COL_ACCESSION)
    sets_acc = [set(dfs[k][COL_ACCESSION]) for k in keys]

    # Unique protein subsets
    unique_dfs = {}
    for idx, k in enumerate(keys):
        others = set()
        for j, kj in enumerate(keys):
            if j != idx: others |= sets_acc[j]
        unique_dfs[k] = dfs[k][dfs[k][COL_ACCESSION].isin(sets_acc[idx] - others)].copy()

    # Core proteome
    core_accs = sets_acc[0] & sets_acc[1] & sets_acc[2] & sets_acc[3]
    core_df = dfs[keys[0]][dfs[keys[0]][COL_ACCESSION].isin(core_accs)].copy()

    # Per-sample enrichment
    enrich = {}
    for k in keys:
        enrich[k] = {}
        for db_name, mem in all_dbs.items():
            enrich[k][db_name] = pathway_fisher_enrichment(dfs[k], bg, mem)

    # Unique protein enrichment
    uniq_enrich = {}
    for k in keys:
        uniq_enrich[k] = {}
        for db_name, mem in all_dbs.items():
            uniq_enrich[k][db_name] = pathway_fisher_enrichment(unique_dfs[k], bg, mem, fdr_thr=0.9) \
                if len(unique_dfs[k]) >= 3 else pd.DataFrame()

    # Core enrichment
    core_enrich = {}
    for db_name, mem in all_dbs.items():
        core_enrich[db_name] = pathway_fisher_enrichment(core_df, bg, mem)

    print(f"  Background: {len(bg)} proteins | Core: {len(core_df)} proteins")
    for k in keys:
        print(f"  {k} unique: {len(unique_dfs[k])} | Reactome: {len(enrich[k]['REACTOME'])} | "
              f"KEGG: {len(enrich[k]['KEGG'])} | WP: {len(enrich[k]['WIKIPATHWAYS'])}")

    return enrich, uniq_enrich, core_enrich, unique_dfs, core_df, bg

# ─────────────────────────────────────────────────────────────────────────────
#  FIGURES 1–19  (descriptive / comparative)
# ─────────────────────────────────────────────────────────────────────────────

def fig1(dfs, meta, keys, shorts, colors):
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle("Figure 1: Individual Sample Protein Identification Summaries\n"
                 "Coral Tissue Proteomics – Gulf of Mannar, Tamil Nadu",
                 fontsize=16, fontweight="bold", y=0.98)
    for idx, k in enumerate(keys):
        ax = axes[idx//2][idx%2]; df = dfs[k]; col = colors[idx]; sp = shorts[idx]
        sc = ax.scatter(df[COL_MW], df[COL_COVERAGE], c=df[COL_PI], cmap="RdYlBu_r",
                        s=df[COL_PEPTIDES]*4+20, alpha=0.75, edgecolors="white", linewidth=0.5)
        cb = plt.colorbar(sc, ax=ax, pad=0.01); cb.set_label("Calc. pI", fontsize=9)
        ax.set_xlabel("Molecular Weight (kDa)", fontsize=10)
        ax.set_ylabel("Sequence Coverage (%)", fontsize=10)
        ax.set_title(f"Sample {idx+1}: $\\it{{{sp}}}$\n"
                     f"n={len(df)} | Avg Cov={df[COL_COVERAGE].mean():.1f}% | "
                     f"MW {df[COL_MW].min():.0f}–{df[COL_MW].max():.0f} kDa",
                     fontsize=10, fontweight="bold", color=col)
        for sz, lbl in [(20,"2 pep"),(60,"10 pep"),(100,"20 pep")]:
            ax.scatter([], [], s=sz, c="gray", alpha=0.6, label=lbl)
        ax.legend(title="# Peptides", fontsize=8, title_fontsize=8, loc="upper right")
        ax.set_facecolor("#f8f9fa"); ax.grid(alpha=0.3)
    plt.tight_layout(rect=[0,0,1,0.96])
    savefig("Fig01_PerSample_Summary.png", fig)
    csv_rows=[]
    for idx,k in enumerate(keys):
        df=dfs[k]
        for _,r in df.iterrows():
            csv_rows.append({"sample":k,"species":meta[k]["species"],
                              "accession":r.get(COL_ACCESSION),"MW_kDa":r.get(COL_MW),
                              "coverage_pct":r.get(COL_COVERAGE),"calc_pI":r.get(COL_PI),
                              "peptides":r.get(COL_PEPTIDES)})
    save_csv("Fig01_PerSample_Summary.csv", pd.DataFrame(csv_rows))

def fig2(dfs, meta, keys, shorts, colors):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Figure 2: Comparative Protein Identification Metrics", fontsize=14, fontweight="bold")
    counts = [len(dfs[k]) for k in keys]
    bars = axes[0].bar(shorts, counts, color=colors, edgecolor="white", linewidth=1.2, width=0.55)
    for bar, val in zip(bars, counts):
        axes[0].text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5, str(val),
                     ha="center", va="bottom", fontweight="bold", fontsize=11)
    axes[0].set_ylabel("Total Proteins"); axes[0].set_title("A. Protein Count", fontweight="bold")
    axes[0].set_xticklabels([f"$\\it{{{s}}}$" for s in shorts], rotation=25, ha="right")
    axes[0].set_facecolor("#f8f9fa"); axes[0].grid(axis="y", alpha=0.3)
    vp = axes[1].violinplot([dfs[k][COL_COVERAGE].dropna().values for k in keys],
                             positions=range(4), widths=0.6, showmedians=True)
    for pc, col in zip(vp["bodies"], colors): pc.set_facecolor(col); pc.set_alpha(0.75)
    vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(2)
    axes[1].set_xticks(range(4))
    axes[1].set_xticklabels([f"$\\it{{{s}}}$" for s in shorts], rotation=25, ha="right")
    axes[1].set_ylabel("Coverage (%)"); axes[1].set_title("B. Coverage Distribution", fontweight="bold")
    axes[1].set_facecolor("#f8f9fa"); axes[1].grid(axis="y", alpha=0.3)
    bp = axes[2].boxplot([dfs[k][COL_PEPTIDES].dropna().values for k in keys],
                          patch_artist=True, widths=0.5, medianprops={"color":"black","linewidth":2})
    for patch, col in zip(bp["boxes"], colors): patch.set_facecolor(col); patch.set_alpha(0.75)
    axes[2].set_xticklabels([f"$\\it{{{s}}}$" for s in shorts], rotation=25, ha="right")
    axes[2].set_ylabel("# Peptides"); axes[2].set_title("C. Peptide Distribution", fontweight="bold")
    axes[2].set_facecolor("#f8f9fa"); axes[2].grid(axis="y", alpha=0.3)
    plt.tight_layout(rect=[0,0,1,0.93]); savefig("Fig02_Comparative_Metrics.png", fig)
    csv_rows=[{"species":meta[k]["species"],"sample":k,"protein_count":len(dfs[k]),
               "coverage_mean":dfs[k][COL_COVERAGE].mean(),"coverage_median":dfs[k][COL_COVERAGE].median(),
               "coverage_std":dfs[k][COL_COVERAGE].std(),
               "peptides_mean":dfs[k][COL_PEPTIDES].mean(),"peptides_median":dfs[k][COL_PEPTIDES].median(),
               "peptides_std":dfs[k][COL_PEPTIDES].std()} for k in keys]
    save_csv("Fig02_Comparative_Metrics.csv", pd.DataFrame(csv_rows))

def fig3(dfs, meta, keys, shorts, colors):
    sets = [set(dfs[k][COL_ACCESSION]) for k in keys]
    all_acc = set().union(*sets)
    combos = []
    for mask in range(1, 16):
        mem=[i for i in range(4) if mask&(1<<i)]; non=[i for i in range(4) if not(mask&(1<<i))]
        present=all_acc.copy()
        for i in mem: present&=sets[i]
        for i in non: present-=sets[i]
        if present: combos.append((mem,len(present)))
    combos.sort(key=lambda x:-x[1]); n=len(combos)
    barcols=["#C0392B" if len(m)==4 else colors[m[0]] if len(m)==1 else "#7F8C8D" for m,_ in combos]
    fig=plt.figure(figsize=(max(20,n*1.1),10))
    fig.suptitle("Figure 3: UpSet Plot – Shared and Unique Proteins Across Four Coral Species",
                 fontsize=14,fontweight="bold",y=0.99)
    ax_bar=fig.add_axes([0.22,0.42,0.74,0.50]); ax_dot=fig.add_axes([0.22,0.10,0.74,0.30])
    ax_tot=fig.add_axes([0.02,0.10,0.18,0.30])
    xs=np.arange(n); bars=ax_bar.bar(xs,[c[1] for c in combos],color=barcols,edgecolor="white",linewidth=0.8)
    for bar,(_,cnt) in zip(bars,combos):
        ax_bar.text(bar.get_x()+bar.get_width()/2,bar.get_height()+0.3,str(cnt),
                    ha="center",va="bottom",fontsize=8,fontweight="bold")
    ax_bar.set_xlim(-0.5,n-0.5); ax_bar.set_ylabel("Protein Count",fontsize=11)
    ax_bar.set_xticks([]); ax_bar.set_facecolor("#f8f9fa"); ax_bar.grid(axis="y",alpha=0.3)
    ax_bar.spines["bottom"].set_visible(False)
    ax_bar.legend(handles=[mpatches.Patch(color="#C0392B",label="Core (all 4)"),
        mpatches.Patch(color=colors[0],label=f"$\\it{{{shorts[0]}}}$ only"),
        mpatches.Patch(color=colors[1],label=f"$\\it{{{shorts[1]}}}$ only"),
        mpatches.Patch(color=colors[2],label=f"$\\it{{{shorts[2]}}}$ only"),
        mpatches.Patch(color=colors[3],label=f"$\\it{{{shorts[3]}}}$ only"),
        mpatches.Patch(color="#7F8C8D",label="Pairwise/Triple")],fontsize=8,loc="upper right",framealpha=0.9)
    ax_dot.set_facecolor("#f8f9fa"); ax_dot.set_xlim(-0.5,n-0.5); ax_dot.set_ylim(-0.5,3.5)
    ax_dot.set_yticks(range(4)); ax_dot.set_yticklabels([f"$\\it{{{s}}}$" for s in reversed(shorts)],fontsize=9)
    ax_dot.set_xticks([])
    for xi,(mem,_) in enumerate(combos):
        for yi in range(4): ax_dot.scatter(xi,3-yi,s=120,color="#d0d0d0",zorder=2)
        pys=[]
        for mi in mem: ax_dot.scatter(xi,3-mi,s=150,color=colors[mi],zorder=3); pys.append(3-mi)
        if len(pys)>1: ax_dot.plot([xi,xi],[min(pys),max(pys)],color="#555",lw=2.5,zorder=1)
    totals=[len(s) for s in sets]
    ax_tot.barh([3,2,1,0],totals,color=colors,edgecolor="white",height=0.5)
    for i,(t,col) in enumerate(zip(totals,colors)):
        ax_tot.text(t+0.5,3-i,str(t),va="center",fontsize=9,fontweight="bold",color=col)
    ax_tot.set_yticks([3,2,1,0]); ax_tot.set_yticklabels([f"$\\it{{{s}}}$" for s in shorts],fontsize=9)
    ax_tot.set_xlabel("Total Proteins",fontsize=10); ax_tot.set_title("Set Size",fontsize=10,fontweight="bold")
    ax_tot.set_facecolor("#f8f9fa"); ax_tot.invert_xaxis(); ax_tot.grid(axis="x",alpha=0.3)
    savefig("Fig03_UpSet_Plot.png",fig)
    csv_rows=[{"species_combo":" + ".join(shorts[i] for i in mem),"n_species_in_combo":len(mem),
               "protein_count":cnt} for mem,cnt in combos]
    csv_rows+=[{"species_combo":f"{shorts[i]} (set total)","n_species_in_combo":1,
                "protein_count":len(sets[i])} for i in range(4)]
    save_csv("Fig03_UpSet_Plot.csv", pd.DataFrame(csv_rows))

def fig4(dfs, meta, keys, shorts, colors):
    mat=np.zeros((4,4))
    for i,ki in enumerate(keys):
        for j,kj in enumerate(keys):
            a=set(dfs[ki][COL_ACCESSION]); b=set(dfs[kj][COL_ACCESSION])
            mat[i,j]=round(len(a&b)/len(a|b),3) if(a|b) else 0
    fig,ax=plt.subplots(figsize=(8,6.5))
    im=ax.imshow(mat,cmap="Blues",vmin=0,vmax=1); plt.colorbar(im,ax=ax,label="Jaccard Similarity",shrink=0.8)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],rotation=35,ha="right",fontsize=10)
    ax.set_yticklabels([f"$\\it{{{s}}}$" for s in shorts],fontsize=10)
    for i in range(4):
        for j in range(4):
            ax.text(j,i,f"{mat[i,j]:.3f}",ha="center",va="center",
                    color="white" if mat[i,j]>0.55 else "black",fontsize=12,fontweight="bold")
    ax.set_title("Figure 4: Pairwise Jaccard Similarity Matrix",fontsize=12,fontweight="bold",pad=12)
    plt.tight_layout(); savefig("Fig04_Jaccard_Heatmap.png",fig)
    jac_df=pd.DataFrame(mat,index=shorts,columns=shorts); jac_df.insert(0,"species",shorts)
    save_csv("Fig04_Jaccard_Heatmap.csv", jac_df)

def fig5(dfs, meta, keys, shorts, colors):
    all_bp=set(); bp_counts={}
    for k in keys: c=parse_go(dfs[k][COL_BP]); bp_counts[k]=c; all_bp|=set(c.keys())
    pooled=Counter()
    for c in bp_counts.values(): pooled.update(c)
    keep, other_key = top_terms_with_other(pooled, TOP_N_GO_PLOT, "Other BP terms")
    plot_terms = keep + ([other_key] if other_key else [])
    palette=plt.cm.tab20(np.linspace(0,1,len(plot_terms)))
    fig,ax=plt.subplots(figsize=(14,8)); bottoms=np.zeros(4)
    for bi,term in enumerate(plot_terms):
        if term==other_key:
            vals=[sum(cnt for t,cnt in bp_counts[k].items() if t not in keep) for k in keys]
        else:
            vals=[bp_counts[k].get(term,0) for k in keys]
        ax.bar(range(4),vals,bottom=bottoms,color=palette[bi],label=term,edgecolor="white",linewidth=0.5)
        bottoms+=np.array(vals)
    ax.set_xticks(range(4)); ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],fontsize=11)
    ax.set_ylabel("Protein Count",fontsize=12)
    ax.set_title(f"Figure 5: GO Biological Process Classification – Stacked Bar\n"
                 f"(top {TOP_N_GO_PLOT} terms by pooled frequency; rest grouped as 'Other')",
                 fontsize=12,fontweight="bold")
    ax.legend(title="Biological Process",bbox_to_anchor=(1.01,1),loc="upper left",fontsize=8,ncol=1)
    ax.set_facecolor("#f8f9fa"); ax.grid(axis="y",alpha=0.3); plt.tight_layout(rect=[0,0,0.78,0.92])
    savefig("Fig05_GO_BP_Stacked.png",fig)
    csv_rows=[{"biological_process":term,**{shorts[i]:bp_counts[k].get(term,0) for i,k in enumerate(keys)}}
              for term in sorted(all_bp)]
    save_csv("Fig05_GO_BP_Stacked.csv", pd.DataFrame(csv_rows))

def fig6(dfs, meta, keys, shorts, colors):
    fig,axes=plt.subplots(1,2,figsize=(20,9))
    fig.suptitle(f"Figure 6: GO Cellular Component & Molecular Function – Bubble Charts\n"
                 f"(top {TOP_N_GO_PLOT} terms per ontology by pooled frequency)",fontsize=13,fontweight="bold")
    for ai,(go_col,label) in enumerate([(COL_CC,"Cellular Component"),(COL_MF,"Molecular Function")]):
        ax=axes[ai]; all_counts={}
        pooled=Counter()
        for k in keys: c=parse_go(dfs[k][go_col]); all_counts[k]=c; pooled.update(c)
        keep,_ = top_terms_with_other(pooled, TOP_N_GO_PLOT)
        all_terms=sorted(keep, key=lambda t:-pooled[t]); y_pos={t:i for i,t in enumerate(all_terms)}
        for xi,(k,col) in enumerate(zip(keys,colors)):
            for term,cnt in all_counts[k].items():
                if term not in y_pos: continue
                ax.scatter(xi,y_pos[term],s=cnt*40+30,color=col,alpha=0.75,edgecolors="white",linewidth=0.8)
                ax.text(xi+0.08,y_pos[term],str(cnt),va="center",fontsize=7,color=col)
        ax.set_xticks(range(4)); ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],rotation=25,ha="right",fontsize=9)
        ax.set_yticks(list(y_pos.values())); ax.set_yticklabels(list(all_terms),fontsize=8)
        ax.set_ylim(-1,len(all_terms))
        ax.set_title(f"GO {label}",fontsize=11,fontweight="bold"); ax.set_xlim(-0.5,3.5)
        ax.set_facecolor("#f8f9fa"); ax.grid(alpha=0.2)
        ax.legend(handles=[mpatches.Patch(color=colors[i],label=f"$\\it{{{shorts[i]}}}$") for i in range(4)],fontsize=8,loc="lower right")
    plt.tight_layout(rect=[0,0,1,0.90]); savefig("Fig06_GO_CC_MF_Bubble.png",fig)
    csv_rows=[]
    for go_col,label in [(COL_CC,"Cellular Component"),(COL_MF,"Molecular Function")]:
        all_counts={k:parse_go(dfs[k][go_col]) for k in keys}
        all_terms=sorted(set().union(*[set(c.keys()) for c in all_counts.values()]))
        for term in all_terms:
            csv_rows.append({"ontology":label,"term":term,
                              **{shorts[i]:all_counts[k].get(term,0) for i,k in enumerate(keys)}})
    save_csv("Fig06_GO_CC_MF_Bubble.csv", pd.DataFrame(csv_rows))

def fig7(dfs, meta, keys, shorts, colors):
    fig,axes=plt.subplots(1,3,figsize=(22,9))
    fig.suptitle("Figure 7: Comparative GO Annotation Heatmaps",fontsize=14,fontweight="bold")
    csv_rows=[]
    for ai,(go_col,cmap,label) in enumerate([(COL_BP,"Blues","Biological Process"),
                                              (COL_CC,"Oranges","Cellular Component"),
                                              (COL_MF,"Greens","Molecular Function")]):
        ax=axes[ai]; all_terms=set(); all_counts={}
        for k in keys: c=parse_go(dfs[k][go_col]); all_counts[k]=c; all_terms|=set(c.keys())
        pooled=Counter()
        for c in all_counts.values(): pooled.update(c)
        keep,_ = top_terms_with_other(pooled, TOP_N_GO_PLOT)
        all_terms=sorted(keep, key=lambda t:-pooled[t])
        if not all_terms:
            ax.text(0.5,0.5,"No terms retrieved\nfor this ontology",ha="center",va="center",
                    fontsize=11,color="gray",transform=ax.transAxes)
            ax.set_title(label,fontsize=10,fontweight="bold"); ax.axis("off")
            continue
        mat=np.array([[all_counts[k].get(t,0) for k in keys] for t in all_terms])
        rs=mat.sum(axis=1,keepdims=True); mat_norm=np.where(rs>0,mat/rs,0)
        im=ax.imshow(mat_norm,cmap=cmap,aspect="auto",vmin=0,vmax=1)
        plt.colorbar(im,ax=ax,label="Proportion",shrink=0.6)
        ax.set_xticks(range(4)); ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],rotation=35,ha="right",fontsize=9)
        ax.set_yticks(range(len(all_terms))); ax.set_yticklabels(all_terms,fontsize=8)
        for i in range(len(all_terms)):
            for j in range(4):
                v=mat[i,j]
                if v>0: ax.text(j,i,str(v),ha="center",va="center",fontsize=7,
                                color="white" if mat_norm[i,j]>0.7 else "black")
        ax.set_title(f"{label}\n(top {TOP_N_GO_PLOT} terms)",fontsize=10,fontweight="bold")
        for term in all_terms:
            csv_rows.append({"ontology":label,"term":term,
                              **{shorts[i]:all_counts[k].get(term,0) for i,k in enumerate(keys)}})
    plt.tight_layout(rect=[0,0,1,0.90]); savefig("Fig07_GO_Comparative_Heatmap.png",fig)
    save_csv("Fig07_GO_Comparative_Heatmap.csv", pd.DataFrame(csv_rows))

def fig8(dfs, meta, keys, shorts, colors):
    flows={}; bp_all=set(); pooled=Counter()
    for si,k in enumerate(keys):
        for _,row in dfs[k].iterrows():
            for bp in [t.strip() for t in str(row.get(COL_BP,"")).split(";") if t.strip() and t.strip()!="nan"]:
                bp_all.add(bp); flows[(si,bp)]=flows.get((si,bp),0)+1; pooled[bp]+=1
    keep, other_key = top_terms_with_other(pooled, TOP_N_GO_PLOT, "Other BP terms")
    plot_flows = dict(flows)
    if other_key:
        for si in range(4):
            other_cnt = sum(cnt for (s2,bp),cnt in flows.items() if s2==si and bp not in keep)
            if other_cnt>0: plot_flows[(si,other_key)]=other_cnt
        for k2 in list(plot_flows.keys()):
            if k2[1] not in keep and k2[1]!=other_key: del plot_flows[k2]
    bp_list = keep + ([other_key] if other_key else [])
    sp_h=[len(dfs[k]) for k in keys]
    total_sp=sum(sp_h); total_bp=sum(plot_flows.values()); gap=0.02; frac=0.85
    bp_cmap=plt.cm.tab20(np.linspace(0,1,max(len(bp_list),1)))
    fig,ax=plt.subplots(figsize=(18,12)); ax.set_facecolor("#f0f0f0"); ax.axis("off")
    ax.set_title(f"Figure 8: Sankey / Alluvial – Species → Biological Process\n"
                 f"(top {TOP_N_GO_PLOT} BP terms by pooled frequency; rest grouped as 'Other')",
                 fontsize=12,fontweight="bold",pad=15)
    x_sp=0.05; x_bp=0.50; sp_y={}; bp_y={}; offsets=[{} for _ in range(4)]
    y=0.07
    for si,k in enumerate(keys):
        h=sp_h[si]/total_sp*frac; sp_y[si]=(y,y+h)
        ax.add_patch(plt.Rectangle((x_sp-0.04,y),0.08,h,color=colors[si],alpha=0.9,zorder=3))
        ax.text(x_sp,y+h/2,f"$\\it{{{shorts[si]}}}$\n({sp_h[si]})",
                ha="center",va="center",fontsize=8.5,fontweight="bold",color="white",zorder=4)
        offsets[si]={"cursor":y,"total":sum(plot_flows.get((si,bp),0) for bp in bp_list),"extent":h}; y+=h+gap
    bp_totals={bp:sum(plot_flows.get((si,bp),0) for si in range(4)) for bp in bp_list}
    y=0.07
    for bi,bp in enumerate(bp_list):
        h=bp_totals[bp]/total_bp*frac if total_bp>0 else 0; bp_y[bi]=(y,y+h)
        ax.add_patch(plt.Rectangle((x_bp-0.04,y),0.08,h,color=bp_cmap[bi],alpha=0.85,zorder=3))
        ax.text(x_bp+0.05,y+h/2,bp,ha="left",va="center",fontsize=8,color="#333",zorder=4); y+=h+gap*0.3
    for si in range(4):
        if offsets[si]["total"]==0: continue
        y_off=offsets[si]["cursor"]
        for bi,bp in enumerate(bp_list):
            cnt=plot_flows.get((si,bp),0)
            if cnt==0: continue
            fs=cnt/total_sp*frac; fb=cnt/total_bp*frac if total_bp>0 else 0
            y2s=bp_y[bi][0]+(bp_totals[bp]-cnt)/total_bp*frac if total_bp>0 else 0
            y2e=bp_y[bi][0]+bp_totals[bp]/total_bp*frac if total_bp>0 else 0
            xs=np.linspace(x_sp+0.04,x_bp-0.04,50)
            ax.fill_between(xs,np.interp(xs,[x_sp+0.04,x_bp-0.04],[y_off,y2s]),
                            np.interp(xs,[x_sp+0.04,x_bp-0.04],[y_off+fs,y2e]),
                            alpha=0.25,color=colors[si],zorder=1); y_off+=fs
    ax.text(x_sp,0.96,"Species",ha="center",fontsize=11,fontweight="bold",color="#222")
    ax.text(x_bp,0.96,"Biological Process",ha="center",fontsize=11,fontweight="bold",color="#222")
    ax.set_xlim(0,0.95); ax.set_ylim(0,1); plt.tight_layout(); savefig("Fig08_Sankey_Species_BP.png",fig)
    csv_rows=[{"species":shorts[si],"biological_process":bp,"protein_count":cnt}
              for (si,bp),cnt in flows.items()]
    save_csv("Fig08_Sankey_Species_BP.csv", pd.DataFrame(csv_rows))

def fig9(dfs, meta, keys, shorts, colors):
    sets=[set(dfs[k][COL_ACCESSION]) for k in keys]
    mat=np.array([[len(sets[i]&sets[j]) for j in range(4)] for i in range(4)])
    fig,ax=plt.subplots(figsize=(11,11),subplot_kw={"aspect":"equal"})
    ax.set_facecolor("#0a0a1a"); fig.patch.set_facecolor("#0a0a1a")
    ax.set_xlim(-1.5,1.5); ax.set_ylim(-1.5,1.5); ax.axis("off")
    ax.set_title("Figure 9: Chord Diagram – Protein Sharing Among Four Coral Species",
                 fontsize=13,fontweight="bold",color="white",pad=15)
    arc_r=1.0; gap_deg=8; total_arc=360-4*gap_deg
    row_sums=[sum(mat[i]) for i in range(4)]; total_w=sum(row_sums)
    angles={}; cum=0
    for i in range(4):
        ext=row_sums[i]/total_w*total_arc if total_w>0 else total_arc/4
        angles[i]=(cum+i*gap_deg,cum+i*gap_deg+ext); cum+=ext
    def d2r(d): return d*math.pi/180
    def arc(r,a1,a2,color,lw=16):
        th=np.linspace(d2r(a1),d2r(a2),100)
        ax.plot(r*np.cos(th),r*np.sin(th),color=color,lw=lw,alpha=0.95,solid_capstyle="round",zorder=5)
    def chord(r,a1s,a1e,a2s,a2e,color,alpha=0.28):
        t=np.linspace(0,1,100); O=np.array([0.0,0.0])
        p0=np.array([r*math.cos(d2r(a1s)),r*math.sin(d2r(a1s))]); p1=np.array([r*math.cos(d2r(a1e)),r*math.sin(d2r(a1e))])
        p2=np.array([r*math.cos(d2r(a2s)),r*math.sin(d2r(a2s))]); p3=np.array([r*math.cos(d2r(a2e)),r*math.sin(d2r(a2e))])
        xs1=(1-t)**2*p0[0]+2*(1-t)*t*O[0]+t**2*p3[0]; ys1=(1-t)**2*p0[1]+2*(1-t)*t*O[1]+t**2*p3[1]
        xs2=(1-t)**2*p1[0]+2*(1-t)*t*O[0]+t**2*p2[0]; ys2=(1-t)**2*p1[1]+2*(1-t)*t*O[1]+t**2*p2[1]
        ax.fill(np.concatenate([xs1,xs2[::-1]]),np.concatenate([ys1,ys2[::-1]]),color=color,alpha=alpha,zorder=1)
    aoff=[{"cursor":angles[i][0],"total":row_sums[i],"extent":angles[i][1]-angles[i][0]} for i in range(4)]
    for i in range(4):
        arc(arc_r,angles[i][0],angles[i][1],colors[i])
        mid=(angles[i][0]+angles[i][1])/2
        ax.text(1.25*math.cos(d2r(mid)),1.25*math.sin(d2r(mid)),
                f"$\\it{{{shorts[i]}}}$\nn={len(sets[i])}",
                ha="center",va="center",fontsize=9.5,color=colors[i],fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3",facecolor="#111",edgecolor=colors[i],alpha=0.9))
    for i in range(4):
        for j in range(i+1,4):
            shared=mat[i,j]
            if shared==0: continue
            ei=shared/aoff[i]["total"]*aoff[i]["extent"] if aoff[i]["total"]>0 else 0
            ej=shared/aoff[j]["total"]*aoff[j]["extent"] if aoff[j]["total"]>0 else 0
            ais=aoff[i]["cursor"]; aie=ais+ei; ajs=aoff[j]["cursor"]; aje=ajs+ej
            aoff[i]["cursor"]+=ei; aoff[j]["cursor"]+=ej
            ci=np.array(mcolors.to_rgba(colors[i])); cj=np.array(mcolors.to_rgba(colors[j]))
            chord(arc_r,ais,aie,ajs,aje,color=tuple((ci+cj)/2))
    ax.text(0,0,"Coral\nProteomics",ha="center",va="center",fontsize=11,color="white",fontweight="bold",alpha=0.6)
    plt.tight_layout(); savefig("Fig09_Chord_Diagram.png",fig)
    mat_df=pd.DataFrame(mat,index=shorts,columns=shorts); mat_df.insert(0,"species",shorts)
    save_csv("Fig09_Chord_Diagram.csv", mat_df)

def fig10(dfs, meta, keys, shorts, colors):
    edges={}; term_totals=Counter()
    for si,k in enumerate(keys):
        for val in dfs[k][COL_BP].dropna():
            for t in str(val).split(";"):
                t=t.strip()
                if t and t!="nan": ek=(shorts[si],t); edges[ek]=edges.get(ek,0)+1; term_totals[t]+=1
    top_terms=[t for t,_ in term_totals.most_common(14)]
    bp_cmap=plt.cm.Set3(np.linspace(0,1,len(top_terms)))
    fig,ax=plt.subplots(figsize=(16,14)); ax.set_facecolor("#f5f5f5")
    ax.set_title("Figure 10: Functional Network – Species × Biological Process",fontsize=13,fontweight="bold",pad=12)
    pos={}
    sp_ys=np.linspace(3.2,-3.2,len(shorts))
    for si,sp in enumerate(shorts): pos[sp]=(-2,sp_ys[si])
    term_ys=np.linspace(5.2,-5.2,len(top_terms))
    for ti,t in enumerate(top_terms): pos[t]=(2.5,term_ys[ti])
    max_w=max((edges.get((sp,t),0) for sp in shorts for t in top_terms),default=1)
    for sp,col in zip(shorts,colors):
        for t in top_terms:
            w=edges.get((sp,t),0)
            if w>0: ax.plot([pos[sp][0],pos[t][0]],[pos[sp][1],pos[t][1]],
                            color=col,alpha=0.35,lw=w/max_w*5+0.5,zorder=1)
    for sp,col,k in zip(shorts,colors,keys):
        ax.scatter(*pos[sp],s=2500,color=col,zorder=4,edgecolors="white",linewidths=2)
        ax.text(pos[sp][0]-0.15,pos[sp][1],f"$\\it{{{sp}}}$\n(n={len(dfs[k])})",
                ha="right",va="center",fontsize=9.5,fontweight="bold",color=col)
    for bi,t in enumerate(top_terms):
        ax.scatter(*pos[t],s=min(400+term_totals[t]*30,2000),color=bp_cmap[bi],
                   zorder=4,edgecolors="#555",linewidths=1)
        ax.text(pos[t][0]+0.1,pos[t][1],t,ha="left",va="center",fontsize=9,color="#333",
                bbox=dict(boxstyle="round,pad=0.2",facecolor="white",alpha=0.7,edgecolor="none"))
    ax.set_xlim(-4.5,6.5); ax.set_ylim(-6,6); ax.axis("off")
    for wv,lbl in [(1,"1 protein"),(5,"5 proteins"),(10,"10+")]:
        ax.plot([],[],color="gray",lw=wv/max_w*5+0.5,label=lbl)
    ax.legend(title="Edge weight",fontsize=9,loc="lower left")
    plt.tight_layout(); savefig("Fig10_Functional_Network.png",fig)
    csv_rows=[{"species":sp,"biological_process":t,"protein_count":edges.get((sp,t),0)}
              for sp in shorts for t in top_terms if edges.get((sp,t),0)>0]
    save_csv("Fig10_Functional_Network.csv", pd.DataFrame(csv_rows))

def fig11(dfs, meta, keys, shorts, colors):
    fig,ax=plt.subplots(figsize=(12,8))
    for si,(k,col) in enumerate(zip(keys,colors)):
        df=dfs[k]
        ax.scatter(df[COL_PI],df[COL_MW],s=df[COL_UNIQ_PEP]*15+20,color=col,alpha=0.65,
                   edgecolors="white",linewidth=0.7,label=f"$\\it{{{shorts[si]}}}$ (n={len(df)})")
    ax.set_xlabel("Calculated pI",fontsize=12); ax.set_ylabel("Molecular Weight (kDa)",fontsize=12)
    ax.set_yscale("log")
    ax.set_title("Figure 11: pI vs MW – All Species\n(bubble = # unique peptides)",fontsize=12,fontweight="bold")
    ax.legend(fontsize=10,framealpha=0.9); ax.axvline(x=7,color="gray",linestyle="--",alpha=0.5)
    ax.set_facecolor("#f8f9fa"); ax.grid(alpha=0.3); plt.tight_layout()
    savefig("Fig11_pI_MW_Scatter.png",fig)
    csv_rows=[]
    for si,k in enumerate(keys):
        for _,r in dfs[k].iterrows():
            csv_rows.append({"species":shorts[si],"accession":r.get(COL_ACCESSION),
                              "calc_pI":r.get(COL_PI),"MW_kDa":r.get(COL_MW),
                              "unique_peptides":r.get(COL_UNIQ_PEP)})
    save_csv("Fig11_pI_MW_Scatter.csv", pd.DataFrame(csv_rows))

def fig12(dfs, meta, keys, shorts, colors):
    fig,axes=plt.subplots(1,2,figsize=(16,6))
    fig.suptitle("Figure 12: Protein Size and pI Distribution Comparisons",fontsize=13,fontweight="bold")
    bins_mw=np.linspace(0,max(dfs[k][COL_MW].max() for k in keys)*1.1,30)
    for si,(k,col) in enumerate(zip(keys,colors)):
        axes[0].hist(dfs[k][COL_MW].dropna(),bins=bins_mw,color=col,alpha=0.5,
                     label=f"$\\it{{{shorts[si]}}}$",edgecolor="white",linewidth=0.3)
        axes[0].axvline(dfs[k][COL_MW].median(),color=col,lw=2,linestyle="--",alpha=0.9)
    axes[0].set_xlabel("Molecular Weight (kDa)",fontsize=11); axes[0].set_ylabel("Count",fontsize=11)
    axes[0].set_title("A. MW Distribution",fontweight="bold"); axes[0].legend(fontsize=9)
    axes[0].set_facecolor("#f8f9fa"); axes[0].grid(alpha=0.3)
    bins_pi=np.linspace(3,12,30)
    for si,(k,col) in enumerate(zip(keys,colors)):
        axes[1].hist(dfs[k][COL_PI].dropna(),bins=bins_pi,color=col,alpha=0.5,
                     label=f"$\\it{{{shorts[si]}}}$",edgecolor="white",linewidth=0.3)
        axes[1].axvline(dfs[k][COL_PI].median(),color=col,lw=2,linestyle="--",alpha=0.9)
    axes[1].axvline(x=7,color="black",lw=1.5,alpha=0.4)
    axes[1].set_xlabel("Calculated pI",fontsize=11); axes[1].set_ylabel("Count",fontsize=11)
    axes[1].set_title("B. pI Distribution",fontweight="bold"); axes[1].legend(fontsize=9)
    axes[1].set_facecolor("#f8f9fa"); axes[1].grid(alpha=0.3)
    plt.tight_layout(rect=[0,0,1,0.93]); savefig("Fig12_MW_pI_Distributions.png",fig)
    csv_rows=[{"species":shorts[si],"sample":k,
               "MW_mean_kDa":dfs[k][COL_MW].mean(),"MW_median_kDa":dfs[k][COL_MW].median(),
               "pI_mean":dfs[k][COL_PI].mean(),"pI_median":dfs[k][COL_PI].median()}
              for si,k in enumerate(keys)]
    save_csv("Fig12_MW_pI_Distributions.csv", pd.DataFrame(csv_rows))

def fig13(dfs, meta, keys, shorts, colors, reactome_mem, kegg_mem, wp_mem):
    categories = {
        "GO annotation (BP/CC/MF)": None,   # handled specially below
        "Reactome pathway":  reactome_mem,
        "KEGG pathway":      kegg_mem,
        "WikiPathways":      wp_mem,
    }
    fig,ax=plt.subplots(figsize=(16,7))
    fig.suptitle("Figure 13: Real Database Annotation Coverage per Species\n"
                 "(proteins with ≥1 curated GO term / pathway hit, fetched live by accession)",
                 fontsize=12,fontweight="bold")
    yp=np.arange(len(categories)); bw=0.18; offs=np.linspace(-1.5*bw,1.5*bw,4)
    csv_rows=[]
    for si,(k,col) in enumerate(zip(keys,colors)):
        accs = list(dfs[k][COL_ACCESSION]); n = len(accs)
        counts=[]
        for cat, mem in categories.items():
            if mem is None:
                c = int(((dfs[k][COL_BP].fillna("")!="") | (dfs[k][COL_CC].fillna("")!="") |
                         (dfs[k][COL_MF].fillna("")!="")).sum())
            else:
                c = sum(1 for a in accs if mem.get(a))
            counts.append(c)
        bars=ax.barh(yp+offs[si],counts,height=bw,color=col,alpha=0.82,
                     label=f"$\\it{{{shorts[si]}}}$",edgecolor="white")
        for bar,val in zip(bars,counts):
            pct = 100*val/n if n else 0
            ax.text(bar.get_width()+max(n*0.01,0.5),bar.get_y()+bar.get_height()/2,
                    f"{val} ({pct:.0f}%)",va="center",fontsize=7.5,color=col)
        for cat,val in zip(categories.keys(),counts):
            csv_rows.append({"species":shorts[si],"database":cat,"proteins_annotated":val,
                              "total_proteins":n,"pct_annotated":round(100*val/n,1) if n else 0})
    ax.set_yticks(yp); ax.set_yticklabels(list(categories.keys()),fontsize=10)
    ax.set_xlabel("Proteins with ≥1 real database hit",fontsize=11); ax.legend(fontsize=10,loc="lower right")
    ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3); plt.tight_layout(rect=[0,0,1,0.88])
    savefig("Fig13_Annotation_Database_Coverage.png",fig)
    save_csv("Fig13_Annotation_Database_Coverage.csv", pd.DataFrame(csv_rows))

def fig14(dfs, meta, keys, shorts, colors):
    sets=[set(dfs[k][COL_ACCESSION]) for k in keys]; all_acc=set().union(*sets)
    presence={}; cov_matrix={}; desc_map={}
    for acc in all_acc:
        presence[acc]=0; cov_matrix[acc]=[]
        for k in keys:
            row=dfs[k][dfs[k][COL_ACCESSION]==acc]
            if len(row)>0:
                presence[acc]+=1; cov_matrix[acc].append(row.iloc[0][COL_COVERAGE])
                desc_map[acc]=short_desc(row.iloc[0][COL_DESC])
            else: cov_matrix[acc].append(0)
    shared=sorted([a for a,p in presence.items() if p>=2],key=lambda a:-presence[a])[:40]
    if not shared: print("  [Fig14] No shared proteins."); return
    mat=np.array([cov_matrix[a] for a in shared])
    fig,ax=plt.subplots(figsize=(12,max(10,len(shared)*0.4)))
    im=ax.imshow(mat,cmap="YlOrRd",aspect="auto",vmin=0,vmax=100)
    plt.colorbar(im,ax=ax,label="Coverage (%)",shrink=0.6,pad=0.02)
    ax.set_xticks(range(4)); ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],fontsize=10,rotation=25,ha="right")
    ax.set_yticks(range(len(shared))); ax.set_yticklabels([desc_map.get(a,a) for a in shared],fontsize=7.5)
    for i in range(len(shared)):
        for j in range(4):
            v=mat[i,j]
            if v>0: ax.text(j,i,f"{v:.0f}%",ha="center",va="center",fontsize=6.5,
                            color="white" if v>55 else "#333")
    ax.set_title(f"Figure 14: Cross-Species Coverage Heatmap – Top {len(shared)} Shared Proteins",
                 fontsize=11,fontweight="bold",pad=12)
    plt.tight_layout(); savefig("Fig14_Shared_Protein_Heatmap.png",fig)
    csv_rows=[{"accession":a,"description":desc_map.get(a,a),"n_species_present":presence[a],
               **{shorts[i]:cov_matrix[a][i] for i in range(4)}} for a in shared]
    save_csv("Fig14_Shared_Protein_Heatmap.csv", pd.DataFrame(csv_rows))

def fig15(dfs, meta, keys, shorts, colors):
    sets=[set(dfs[k][COL_ACCESSION]) for k in keys]
    metrics=["Total\nProteins","Avg Coverage\n(%)","Avg Peptides","Unique\nProteins","Avg MW\n(kDa/10)"]
    vals_raw=[]
    for si,k in enumerate(keys):
        df=dfs[k]; others=set()
        for j,kj in enumerate(keys):
            if j!=si: others|=sets[j]
        vals_raw.append([len(df),df[COL_COVERAGE].mean(),df[COL_PEPTIDES].mean(),
                         len(sets[si]-others),df[COL_MW].mean()/10])
    vals_raw=np.array(vals_raw); vmin=vals_raw.min(axis=0); vmax=vals_raw.max(axis=0)
    vals_norm=np.where(vmax-vmin>0,(vals_raw-vmin)/(vmax-vmin),0.5)
    n=len(metrics); angles=np.linspace(0,2*np.pi,n,endpoint=False).tolist(); angles+=angles[:1]
    fig,ax=plt.subplots(figsize=(10,10),subplot_kw={"polar":True}); ax.set_facecolor("#f8f9fa")
    for si in range(4):
        v=vals_norm[si].tolist()+[vals_norm[si][0]]
        ax.plot(angles,v,color=colors[si],lw=2.5,label=f"$\\it{{{shorts[si]}}}$",zorder=3)
        ax.fill(angles,v,color=colors[si],alpha=0.15,zorder=2)
        for ai,(a,raw) in enumerate(zip(angles[:-1],vals_raw[si])):
            ax.annotate(f"{raw:.1f}",xy=(a,vals_norm[si][ai]),fontsize=7.5,color=colors[si],
                        ha="center",va="bottom",xytext=(a,vals_norm[si][ai]+0.08),textcoords="data")
    ax.set_xticks(angles[:-1]); ax.set_xticklabels(metrics,fontsize=10,fontweight="bold")
    ax.set_ylim(0,1.2); ax.set_yticks([0.25,0.5,0.75,1.0])
    ax.set_yticklabels(["25%","50%","75%","100%"],fontsize=7,color="gray")
    ax.set_title("Figure 15: Radar Chart – Normalised Proteomic Profile",fontsize=12,fontweight="bold",pad=25)
    ax.legend(loc="lower right",bbox_to_anchor=(1.35,-0.1),fontsize=10); ax.grid(color="gray",alpha=0.3)
    plt.tight_layout(); savefig("Fig15_Radar_Proteome_Profile.png",fig)
    metric_names=[m.replace("\n"," ") for m in metrics]
    csv_rows=[{"species":shorts[si],**dict(zip(metric_names,vals_raw[si]))} for si in range(4)]
    save_csv("Fig15_Radar_Proteome_Profile.csv", pd.DataFrame(csv_rows))

def fig16(dfs, meta, keys, shorts, colors):
    fig,axes=plt.subplots(2,2,figsize=(18,16))
    fig.suptitle(f"Figure 16: GO Biological Process – Pie Charts per Species\n"
                 f"(top {TOP_N_GO_PLOT} terms per species; rest grouped as 'Other')",
                 fontsize=13,fontweight="bold")
    csv_rows=[]
    for idx,(k,col,sp) in enumerate(zip(keys,colors,shorts)):
        ax=axes[idx//2][idx%2]; counts=parse_go(dfs[k][COL_BP])
        if counts:
            keep, other_key = top_terms_with_other(counts, TOP_N_GO_PLOT, "Other")
            labels = sorted(keep, key=lambda t:-counts[t])
            sizes = [counts[t] for t in labels]
            if other_key:
                labels = labels + [other_key]
                sizes = sizes + [sum(c for t,c in counts.items() if t not in keep)]
            pie_cmap=plt.cm.tab20(np.linspace(0,1,len(labels)))
            wedges,_,ats=ax.pie(sizes,labels=None,autopct="%1.1f%%",colors=pie_cmap,
                                startangle=90,wedgeprops={"edgecolor":"white","linewidth":1.5},pctdistance=0.75)
            for at in ats: at.set_fontsize(8)
            # Legend placed to the side (not below) so it never collides with the row beneath it
            ax.legend(wedges,labels,loc="center left",bbox_to_anchor=(1.02,0.5),ncol=1,fontsize=7.5,framealpha=0.85)
            for lbl,sz in zip(list(counts.keys()),list(counts.values())):
                csv_rows.append({"species":sp,"biological_process":lbl,"protein_count":sz})
        ax.set_title(f"$\\it{{{sp}}}$\n(n={len(dfs[k])})",fontsize=11,fontweight="bold",color=col,pad=10)
    plt.tight_layout(rect=[0,0,1,0.90])
    fig.subplots_adjust(wspace=0.55, hspace=0.4)
    savefig("Fig16_GO_BP_Pie_Charts.png",fig)
    save_csv("Fig16_GO_BP_Pie_Charts.csv", pd.DataFrame(csv_rows))

def fig17(dfs, meta, keys, shorts, colors):
    sets=[set(dfs[k][COL_ACCESSION]) for k in keys]
    fig,axes=plt.subplots(2,2,figsize=(20,16))
    fig.suptitle("Figure 17: Species-Unique Proteins – Top 15 by Coverage",fontsize=13,fontweight="bold")
    csv_rows=[]
    for idx,(k,col,sp) in enumerate(zip(keys,colors,shorts)):
        ax=axes[idx//2][idx%2]; others=set()
        for j,kj in enumerate(keys):
            if j!=idx: others|=sets[j]
        unique=sets[idx]-others; n_unique=len(unique)
        df_uniq=dfs[k][dfs[k][COL_ACCESSION].isin(unique)].copy()
        df_uniq=df_uniq.nlargest(min(15,len(df_uniq)),COL_COVERAGE)
        if len(df_uniq)==0:
            ax.text(0.5,0.5,"No species-unique proteins",ha="center",va="center",transform=ax.transAxes)
            ax.set_title(f"$\\it{{{sp}}}$",color=col,fontweight="bold"); continue
        df_uniq["sd"]=df_uniq[COL_DESC].apply(short_desc)
        bars=ax.barh(range(len(df_uniq)),df_uniq[COL_COVERAGE].values,color=col,alpha=0.80,edgecolor="white")
        ax.set_yticks(range(len(df_uniq))); ax.set_yticklabels(df_uniq["sd"].values,fontsize=8.5)
        for bar,val in zip(bars,df_uniq[COL_COVERAGE].values):
            ax.text(bar.get_width()+0.4,bar.get_y()+bar.get_height()/2,f"{val:.1f}%",va="center",fontsize=8)
        ax.set_xlabel("Coverage (%)",fontsize=10)
        ax.set_title(f"$\\it{{{sp}}}$ – {n_unique} Unique Proteins",fontsize=10,fontweight="bold",color=col)
        ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3)
        for _,r in df_uniq.iterrows():
            csv_rows.append({"species":sp,"accession":r.get(COL_ACCESSION),"description":r.get(COL_DESC),
                              "coverage_pct":r.get(COL_COVERAGE),"total_species_unique_count":n_unique})
    plt.tight_layout(rect=[0,0,1,0.94]); savefig("Fig17_Unique_Proteins.png",fig)
    save_csv("Fig17_Unique_Proteins.csv", pd.DataFrame(csv_rows))

def fig18(dfs, meta, keys, shorts, colors):
    sets=[set(dfs[k][COL_ACCESSION]) for k in keys]
    def uniq(k):
        ki=keys.index(k); o=set()
        [o.__ior__(sets[j]) for j in range(4) if j!=ki]; return len(sets[ki]-o)
    defs={
        "Total Proteins":             lambda k: f"{len(dfs[k])}",
        "Avg Coverage (%)":           lambda k: f"{dfs[k][COL_COVERAGE].mean():.1f} ± {dfs[k][COL_COVERAGE].std():.1f}",
        "Avg # Peptides":             lambda k: f"{dfs[k][COL_PEPTIDES].mean():.1f} ± {dfs[k][COL_PEPTIDES].std():.1f}",
        "Avg # Unique Peptides":      lambda k: f"{dfs[k][COL_UNIQ_PEP].mean():.1f} ± {dfs[k][COL_UNIQ_PEP].std():.1f}",
        "Avg # PSMs":                 lambda k: f"{dfs[k][COL_PSM].mean():.1f} ± {dfs[k][COL_PSM].std():.1f}",
        "Avg MW (kDa)":               lambda k: f"{dfs[k][COL_MW].mean():.1f} ± {dfs[k][COL_MW].std():.1f}",
        "MW Range (kDa)":             lambda k: f"{dfs[k][COL_MW].min():.1f} – {dfs[k][COL_MW].max():.1f}",
        "Avg calc. pI":               lambda k: f"{dfs[k][COL_PI].mean():.2f} ± {dfs[k][COL_PI].std():.2f}",
        "Acidic (pI < 7)":            lambda k: f"{(dfs[k][COL_PI]<7).sum()} ({(dfs[k][COL_PI]<7).mean()*100:.1f}%)",
        "Basic (pI > 7)":             lambda k: f"{(dfs[k][COL_PI]>7).sum()} ({(dfs[k][COL_PI]>7).mean()*100:.1f}%)",
        "Species-unique proteins":    lambda k: f"{uniq(k)}",
        "With GO-BP annotation":      lambda k: f"{dfs[k][COL_BP].notna().sum()} ({dfs[k][COL_BP].notna().mean()*100:.1f}%)",
        "With GO-CC annotation":      lambda k: f"{dfs[k][COL_CC].notna().sum()} ({dfs[k][COL_CC].notna().mean()*100:.1f}%)",
        "With GO-MF annotation":      lambda k: f"{dfs[k][COL_MF].notna().sum()} ({dfs[k][COL_MF].notna().mean()*100:.1f}%)",
    }
    rows=[[metric]+[fn(k) for k in keys] for metric,fn in defs.items()]
    col_labels=["Metric"]+[f"$\\it{{{s}}}$\n(S{i+1})" for i,s in enumerate(shorts)]
    fig,ax=plt.subplots(figsize=(18,8)); ax.axis("off")
    ax.set_title("Figure 18: Comprehensive Proteomic Summary Statistics",fontsize=14,fontweight="bold",pad=15)
    tbl=ax.table(cellText=rows,colLabels=col_labels,cellLoc="center",loc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1,1.55)
    for j in range(5):
        tbl[0,j].set_facecolor("#1B2631"); tbl[0,j].set_text_props(color="white",fontweight="bold")
    bgs=["white","#D6EAF8","#D5F5E3","#FDEBD0","#E8DAEF"]
    for i in range(1,len(rows)+1):
        for j in range(5):
            tbl[i,j].set_facecolor(bgs[j])
            if i%2==0:
                r,g,b,_=mcolors.to_rgba(bgs[j])
                tbl[i,j].set_facecolor((max(r-0.04,0),max(g-0.04,0),max(b-0.04,0),1))
        tbl[i,0].set_text_props(fontweight="bold")
    plt.tight_layout(); savefig("Fig18_Summary_Table.png",fig)
    save_csv("Fig18_Summary_Table.csv", pd.DataFrame(rows, columns=["Metric"]+list(keys)))

def fig19(meta, keys, shorts, colors):
    fig,ax=plt.subplots(figsize=(20,14)); ax.set_xlim(0,20); ax.set_ylim(0,14)
    ax.axis("off"); ax.set_facecolor("#0D1117"); fig.patch.set_facecolor("#0D1117")
    ax.set_title("Figure 19: Comprehensive Proteomic Analysis Workflow\nCoral Tissue LC-MS/MS – Gulf of Mannar, Tamil Nadu",
                 fontsize=15,fontweight="bold",color="white",pad=10)
    def box(x,y,w,h,txt,bc,tc="white",fs=9.5):
        ax.add_patch(mpatches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.1",
                                              facecolor=bc,edgecolor="white",linewidth=1.2,alpha=0.92))
        ax.text(x+w/2,y+h/2,txt,ha="center",va="center",color=tc,fontsize=fs,fontweight="bold")
    def arr(x1,y1,x2,y2):
        ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="->",color="#58A6FF",lw=2.0))
    for x,txt,col in [(0.3,"🪸 Coral Tissue\nCollection\nErwadi, Gulf of Mannar","#1B4F72"),
                       (4.2,"❄ Storage −80°C\nSterile conditions\nNo bio-replicates","#154360"),
                       (8.1,"🧊 Lyophilization\n50 mg powder\nBead-beat","#1A5276"),
                       (12.0,"🧪 Protein Extraction\n8M Urea / AmBicarb\nBCA Quantification","#1F618D"),
                       (15.9,"✂ Tryptic Digestion\nDTT→IAA→Trypsin\nC18 cleanup","#2471A3")]:
        box(x,12.0,3.5,1.2,txt,col)
    for x in [3.8,7.7,11.7,15.5]: arr(x,12.6,x+0.4,12.6)
    box(0.3,10.0,3.5,1.2,"🔬 LC-MS/MS\nSequest HT\nFDR < 1%","#145A32")
    box(4.2,10.0,3.5,1.2,"📊 Protein Inference\nMaster Proteins\nHigh Confidence","#1E8449")
    box(8.1,10.0,3.5,1.2,"🧬 Proteome Discoverer\nProtein+Peptide lists\n4 Samples","#239B56")
    box(12.0,10.0,7.5,1.2,f"S1: {shorts[0]}  |  S2: {shorts[1]}\nS3: {shorts[2]}  |  S4: {shorts[3]}","#27AE60",fs=8.5)
    for x in [3.8,7.7,11.7]: arr(x,10.6,x+0.4,10.6)
    ax.annotate("",xy=(10,10.0),xytext=(10,12.0),arrowprops=dict(arrowstyle="->",color="#58A6FF",lw=2.5))
    ax.text(10.15,11.0,"BRIC-RGCB\nMass Spec Core",color="#58A6FF",fontsize=8.5,va="center")
    y3=7.8
    for x,txt,col in [(0.3,"Per-Sample\nSummary\n(Fig 1,2)","#6C3483"),(3.0,"UpSet Plot\n(Fig 3)","#7D3C98"),
                       (5.7,"Jaccard\n(Fig 4)","#884EA0"),(8.4,"GO Analysis\n(Fig 5-7,16)","#9B59B6"),
                       (11.1,"CC/MF Bubble\n(Fig 6)","#A569BD"),(13.8,"Heatmaps\n(Fig 14)","#BB8FCE"),
                       (16.5,"Radar\n(Fig 15)","#C39BD3")]:
        box(x,y3,2.5,1.5,txt,col)
    y4=5.6
    for x,txt,col in [(0.3,"Sankey\n(Fig 8)","#784212"),(3.0,"Chord\n(Fig 9,30)","#935116"),
                       (5.7,"Network\n(Fig 10)","#A04000"),(8.4,"pI/MW\n(Fig 11,12)","#BA4A00"),
                       (11.1,"Enrichment\n(Fig 20-27)","#D35400"),(13.8,"Core/PPI\n(Fig 28-31)","#E59866"),
                       (16.5,"Tables\n(Fig 18)","#F0B27A")]:
        box(x,y4,2.5,1.5,txt,col)
    xcs=[1.55,4.25,6.95,9.65,12.35,15.05,17.75]
    for xc in xcs: arr(xc,10.0,xc,y3+1.5); arr(xc,y3,xc,y4+1.5)
    box(0.3,3.6,8.8,1.5,"📄 31 Publication-Quality Figures\nIndividual · Comparative · Functional · Enrichment · Network\nPNG 180 dpi","#1C2833")
    box(9.4,3.6,10.2,1.5,"📊 Quantitative Results\nProtein Lists · GO · Pathway Enrichment (R/K/WP)\nJaccard · Core Proteome · PPI Networks","#1C2833")
    ax.plot([xcs[0],xcs[-1]],[5.2,5.2],color="#58A6FF",lw=1.5,linestyle="--",alpha=0.5)
    for xc in xcs: arr(xc,y4,xc,5.2)
    arr(4.75,5.2,4.75,3.6+1.5); arr(14.45,5.2,14.45,3.6+1.5)
    for si,(col,sp) in enumerate(zip(colors,shorts)):
        ax.add_patch(mpatches.Rectangle((0.3+si*3.5,2.5),0.4,0.4,color=col))
        ax.text(0.8+si*3.5,2.7,f"S{si+1}: $\\it{{{sp}}}$",color="white",fontsize=9,va="center")
    ax.text(10.0,0.8,"Mass Spectrometry Core Facility BRIC-RGCB, Thiruvananthapuram  |  VIT Vellore Bioinformatics  |  Gulf of Mannar Coral Study 2026",
            ha="center",va="center",color="#888",fontsize=8.5)
    plt.tight_layout(); savefig("Fig19_Analysis_Workflow.png",fig)

# ─────────────────────────────────────────────────────────────────────────────
#  FIGURES 20–31  (enrichment / PPI / chord)
# ─────────────────────────────────────────────────────────────────────────────

def fig_enrichment_bars(enrich_results, keys, shorts, colors, db_name, fig_num, title_db, dfs):
    fig,axes=plt.subplots(2,2,figsize=(20,17))
    fig.suptitle(f"Figure {fig_num}: {title_db} Pathway Enrichment per Coral Species\n"
                 f"(Fisher's Exact Test, BH-FDR corrected)",fontsize=14,fontweight="bold")
    csv_rows=[]
    for idx,k in enumerate(keys):
        ax=axes[idx//2][idx%2]; col=colors[idx]; sp=shorts[idx]
        res=enrich_results[k][db_name]
        if len(res)==0:
            ax.text(0.5,0.5,"No significant enrichment\n(q-threshold not reached)",
                    ha="center",va="center",transform=ax.transAxes,fontsize=11,color="gray")
            ax.set_title(f"$\\it{{{sp}}}$",color=col,fontweight="bold"); ax.axis("off"); continue
        res=res.head(12); labels=[pw_short(p) for p in res["pathway"]]
        neg_log=-np.log10(res["padj"].clip(1e-10))
        bars=ax.barh(range(len(res)),neg_log,color=col,alpha=0.80,edgecolor="white",linewidth=0.6)
        ax.set_yticks(range(len(res))); ax.set_yticklabels(labels,fontsize=8)
        for bar,hits in zip(bars,res["study_hits"]):
            ax.text(bar.get_width()+0.05,bar.get_y()+bar.get_height()/2,f"n={hits}",va="center",fontsize=7.5,color=col)
        ax.axvline(-np.log10(0.05),color="red",lw=1.5,linestyle="--",alpha=0.7,label="q=0.05")
        ax.set_xlabel("−log₁₀(q-value)",fontsize=10); ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3)
        ax.legend(fontsize=8,loc="lower right")
        ax.set_title(f"$\\it{{{sp}}}$  (n={len(dfs[k])}, {len(res)} pathways)",color=col,fontweight="bold",fontsize=10)
        for _,r in res.iterrows():
            csv_rows.append({"species":sp,"pathway":r["pathway"],"study_hits":r["study_hits"],
                              "bg_hits":r["bg_hits"],"study_n":r["study_n"],"bg_n":r["bg_n"],
                              "pval":r["pval"],"padj":r["padj"],"ratio":r["ratio"]})
    plt.tight_layout(rect=[0,0,1,0.93])
    savefig(f"Fig{fig_num:02d}_{db_name}_Enrichment_PerSpecies.png",fig)
    save_csv(f"Fig{fig_num:02d}_{db_name}_Enrichment_PerSpecies.csv", pd.DataFrame(csv_rows))

def fig_comparative_heatmap(enrich_results, keys, shorts, colors, db_name, fig_num, title_db):
    all_pw=set()
    for k in keys:
        res=enrich_results[k][db_name]
        if len(res)>0: all_pw|=set(res["pathway"])
    all_pw=sorted(all_pw)
    if not all_pw: print(f"  [Fig{fig_num}] No enriched pathways for {db_name}."); return
    mat=np.zeros((len(all_pw),4)); hit_mat=np.zeros((len(all_pw),4),dtype=int)
    for ci,k in enumerate(keys):
        res=enrich_results[k][db_name]
        for ri,pw in enumerate(all_pw):
            row=res[res["pathway"]==pw]
            if len(row)>0:
                mat[ri,ci]=-np.log10(float(row.iloc[0]["padj"])); hit_mat[ri,ci]=int(row.iloc[0]["study_hits"])
    n_pw=len(all_pw)
    fig,ax=plt.subplots(figsize=(12,max(8,n_pw*0.55+2)))
    im=ax.imshow(mat,cmap="YlOrRd",aspect="auto",vmin=0,vmax=max(mat.max(),0.1))
    plt.colorbar(im,ax=ax,label="−log₁₀(q-value)",shrink=0.6,pad=0.02)
    ax.set_xticks(range(4)); ax.set_xticklabels([f"$\\it{{{s}}}$" for s in shorts],rotation=30,ha="right",fontsize=10)
    ax.set_yticks(range(n_pw)); ax.set_yticklabels([pw_short(p) for p in all_pw],fontsize=8)
    for ri in range(n_pw):
        for ci in range(4):
            h=hit_mat[ri,ci]
            if h>0: ax.text(ci,ri,f"{mat[ri,ci]:.1f}\n(n={h})",ha="center",va="center",fontsize=6.5,
                            color="white" if mat[ri,ci]>mat.max()*0.6 else "#333")
    ax.set_title(f"Figure {fig_num}: Comparative {title_db} Enrichment Heatmap\n"
                 f"(values = −log₁₀ q-value; n = protein hits)",fontsize=12,fontweight="bold",pad=15)
    plt.tight_layout(); savefig(f"Fig{fig_num:02d}_Comparative_{db_name}_Heatmap.png",fig)
    csv_rows=[]
    for ri,pw in enumerate(all_pw):
        row={"pathway":pw}
        for ci,k in enumerate(keys):
            row[f"{shorts[ci]}_neglog10_padj"]=mat[ri,ci]
            row[f"{shorts[ci]}_study_hits"]=hit_mat[ri,ci]
        csv_rows.append(row)
    save_csv(f"Fig{fig_num:02d}_Comparative_{db_name}_Heatmap.csv", pd.DataFrame(csv_rows))

def fig26_unique_go(dfs, unique_dfs, keys, shorts, colors):
    fig,axes=plt.subplots(4,3,figsize=(22,28))
    fig.suptitle("Figure 26: Species-Specific GO Enrichment of Unique Proteins\n"
                 "(BP | CC | MF) – Fisher's Exact Test, BH-FDR",fontsize=13,fontweight="bold")
    csv_rows=[]
    go_cols=[(COL_BP,"GO-BP"),(COL_CC,"GO-CC"),(COL_MF,"GO-MF")]
    for ri,(k,col,sp) in enumerate(zip(keys,colors,shorts)):
        for ci,(go_col,go_lbl) in enumerate(go_cols):
            ax=axes[ri][ci]
            res=go_fisher_enrichment(unique_dfs[k],dfs[k],go_col)
            if len(res)==0:
                ax.text(0.5,0.5,f"No significant\n{go_lbl} enrichment",ha="center",va="center",
                        transform=ax.transAxes,fontsize=10,color="gray")
                ax.set_title(f"$\\it{{{sp}}}$ – {go_lbl}",color=col,fontsize=9,fontweight="bold")
                ax.axis("off"); continue
            res=res.head(10); neg_log=-np.log10(res["padj"].clip(1e-10))
            bars=ax.barh(range(len(res)),neg_log,color=col,alpha=0.78,edgecolor="white",linewidth=0.5)
            ax.set_yticks(range(len(res))); ax.set_yticklabels(res["term"].tolist(),fontsize=7.5)
            for bar,hits in zip(bars,res["uniq_hits"]):
                ax.text(bar.get_width()+0.03,bar.get_y()+bar.get_height()/2,f"n={hits}",va="center",fontsize=7,color=col)
            ax.axvline(-np.log10(0.05),color="red",lw=1.2,linestyle="--",alpha=0.6)
            ax.set_xlabel("−log₁₀(q)",fontsize=8); ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3)
            ax.set_title(f"$\\it{{{sp}}}$ – {go_lbl}\n({len(unique_dfs[k])} unique proteins)",
                         color=col,fontsize=9,fontweight="bold")
            for _,r in res.iterrows():
                csv_rows.append({"species":sp,"ontology":go_lbl,"term":r["term"],
                                  "uniq_hits":r["uniq_hits"],"full_hits":r["full_hits"],
                                  "pval":r["pval"],"padj":r["padj"]})
    plt.tight_layout(rect=[0,0,1,0.965])
    fig.subplots_adjust(hspace=0.55, wspace=0.45)
    savefig("Fig26_UniqueProtein_GO_Enrichment.png",fig)
    save_csv("Fig26_UniqueProtein_GO_Enrichment.csv", pd.DataFrame(csv_rows))

def fig27_unique_pathway(uniq_enrich, unique_dfs, keys, shorts, colors):
    fig,axes=plt.subplots(4,3,figsize=(22,28))
    fig.suptitle("Figure 27: Species-Specific Pathway Enrichment of Unique Proteins\n"
                 "(Reactome | KEGG | WikiPathways) – Fisher's Exact, BH-FDR",fontsize=13,fontweight="bold")
    csv_rows=[]
    dbs=[("REACTOME","Reactome"),("KEGG","KEGG"),("WIKIPATHWAYS","WikiPathways")]
    for ri,(k,col,sp) in enumerate(zip(keys,colors,shorts)):
        for ci,(db_name,db_lbl) in enumerate(dbs):
            ax=axes[ri][ci]; res=uniq_enrich[k][db_name]
            if len(res)==0:
                ax.text(0.5,0.5,f"No significant\n{db_lbl} enrichment",ha="center",va="center",
                        transform=ax.transAxes,fontsize=10,color="gray")
                ax.set_title(f"$\\it{{{sp}}}$ – {db_lbl}",color=col,fontsize=9,fontweight="bold")
                ax.axis("off"); continue
            res=res.head(10); labels=[pw_short(p) for p in res["pathway"]]
            neg_log=-np.log10(res["padj"].clip(1e-10))
            bars=ax.barh(range(len(res)),neg_log,color=col,alpha=0.78,edgecolor="white",linewidth=0.5)
            ax.set_yticks(range(len(res))); ax.set_yticklabels(labels,fontsize=7.5)
            for bar,hits in zip(bars,res["study_hits"]):
                ax.text(bar.get_width()+0.03,bar.get_y()+bar.get_height()/2,f"n={hits}",va="center",fontsize=7,color=col)
            ax.axvline(-np.log10(0.05),color="red",lw=1.2,linestyle="--",alpha=0.6)
            ax.set_xlabel("−log₁₀(q)",fontsize=8); ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3)
            ax.set_title(f"$\\it{{{sp}}}$ – {db_lbl}",color=col,fontsize=9,fontweight="bold")
            for _,r in res.iterrows():
                csv_rows.append({"species":sp,"database":db_lbl,"pathway":r["pathway"],
                                  "study_hits":r["study_hits"],"bg_hits":r["bg_hits"],
                                  "pval":r["pval"],"padj":r["padj"],"ratio":r["ratio"]})
    plt.tight_layout(rect=[0,0,1,0.965])
    fig.subplots_adjust(hspace=0.55, wspace=0.45)
    savefig("Fig27_UniqueProtein_Pathway_Enrichment.png",fig)
    save_csv("Fig27_UniqueProtein_Pathway_Enrichment.csv", pd.DataFrame(csv_rows))

def fig28_core_enrichment(core_enrich, core_df):
    fig,axes=plt.subplots(1,3,figsize=(22,9))
    fig.suptitle("Figure 28: Core Proteome (Shared All 4 Species) – GO + Pathway Enrichment",
                 fontsize=13,fontweight="bold")
    csv_go_rows=[]; csv_pw_rows=[]
    bp_c=parse_go(core_df[COL_BP]); cc_c=parse_go(core_df[COL_CC]); mf_c=parse_go(core_df[COL_MF])
    ax=axes[0]; all_go_terms={}
    for t,n in bp_c.items(): all_go_terms[f"[BP] {t}"]=n
    for t,n in cc_c.items(): all_go_terms[f"[CC] {t}"]=n
    for t,n in mf_c.items(): all_go_terms[f"[MF] {t}"]=n
    sorted_go=sorted(all_go_terms.items(),key=lambda x:-x[1])[:20]
    terms=[x[0] for x in sorted_go]; vals=[x[1] for x in sorted_go]
    go_colors=["#1B4F72" if t.startswith("[BP]") else "#1E8449" if t.startswith("[CC]") else "#784212" for t in terms]
    bars=ax.barh(range(len(terms)),vals,color=go_colors,alpha=0.82,edgecolor="white")
    ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms,fontsize=8)
    for bar,val in zip(bars,vals): ax.text(bar.get_width()+0.05,bar.get_y()+bar.get_height()/2,str(val),va="center",fontsize=8)
    ax.set_xlabel("Count",fontsize=10); ax.set_title(f"A. GO Annotation\n(n={len(core_df)} core proteins)",fontweight="bold")
    ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3)
    ax.legend(handles=[mpatches.Patch(color="#1B4F72",label="BP"),mpatches.Patch(color="#1E8449",label="CC"),
                        mpatches.Patch(color="#784212",label="MF")],fontsize=9)
    csv_go_rows=[{"term":t,"count":v} for t,v in sorted_go]
    for ax_idx,(db_name,db_lbl,panel) in enumerate([("REACTOME","Reactome","B"),("KEGG","KEGG","C")]):
        ax=axes[ax_idx+1]; res=core_enrich[db_name]
        if len(res)==0:
            ax.text(0.5,0.5,"No enriched pathways",ha="center",va="center",transform=ax.transAxes); ax.axis("off"); continue
        res=res.head(14); labels=[pw_short(p) for p in res["pathway"]]
        neg_log=-np.log10(res["padj"].clip(1e-10))
        bars=ax.barh(range(len(res)),neg_log,color="#C0392B",alpha=0.82,edgecolor="white",linewidth=0.6)
        ax.set_yticks(range(len(res))); ax.set_yticklabels(labels,fontsize=8)
        for bar,hits in zip(bars,res["study_hits"]):
            ax.text(bar.get_width()+0.03,bar.get_y()+bar.get_height()/2,f"n={hits}",va="center",fontsize=7.5,color="#C0392B")
        ax.axvline(-np.log10(0.05),color="red",lw=1.5,linestyle="--",alpha=0.7,label="q=0.05")
        ax.set_xlabel("−log₁₀(q)",fontsize=10); ax.set_title(f"{panel}. {db_lbl}\n(Core proteome)",fontweight="bold")
        ax.set_facecolor("#f8f9fa"); ax.grid(axis="x",alpha=0.3); ax.legend(fontsize=8)
        for _,r in res.iterrows():
            csv_pw_rows.append({"database":db_lbl,"pathway":r["pathway"],"study_hits":r["study_hits"],
                                 "pval":r["pval"],"padj":r["padj"]})
    plt.tight_layout(rect=[0,0,1,0.90])
    savefig("Fig28_CoreProteome_GO_Pathway_Enrichment.png",fig)
    save_csv("Fig28_CoreProteome_GO_annotation.csv", pd.DataFrame(csv_go_rows))
    save_csv("Fig28_CoreProteome_Pathway_Enrichment.csv", pd.DataFrame(csv_pw_rows))

def fig29_ppi_core(core_df, pw_mem):
    accs=list(core_df[COL_ACCESSION]); names={r[COL_ACCESSION]:pname(r[COL_DESC]) for _,r in core_df.iterrows()}
    pws={acc:pw_mem.get(acc, set()) for acc in accs}
    gos={r[COL_ACCESSION]:get_go_set(r) for _,r in core_df.iterrows()}
    G=nx.Graph()
    for acc in accs: G.add_node(acc,label=names[acc][:28],pw=pws[acc],go=gos.get(acc,set()))
    for i in range(len(accs)):
        for j in range(i+1,len(accs)):
            a,b=accs[i],accs[j]
            w=len(G.nodes[a]["pw"]&G.nodes[b]["pw"])*2+len(G.nodes[a]["go"]&G.nodes[b]["go"])
            if w>0: G.add_edge(a,b,weight=w)
    print(f"  PPI core: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    from collections import Counter
    def dom_bp(acc):
        go=G.nodes[acc]["go"]; bp=[t for t in go if not t.startswith("other")]
        return bp[0] if bp else "other"
    all_dom=[dom_bp(a) for a in accs]; term_counts=Counter(all_dom)
    top_t=[t for t,_ in term_counts.most_common(8)]
    tc=plt.cm.Set1(np.linspace(0,0.9,len(top_t))); tcmap={t:tc[i] for i,t in enumerate(top_t)}
    pos=nx.spring_layout(G,k=2.8,seed=42,weight="weight",iterations=100)
    fig,ax=plt.subplots(figsize=(20,18)); ax.set_facecolor("#0d0d1a"); fig.patch.set_facecolor("#0d0d1a")
    ax.set_title("Figure 29: PPI Co-association Network – Core Coral Proteome\n"
                 "(Nodes = conserved proteins; edges = shared functional annotations)",
                 fontsize=13,fontweight="bold",color="white",pad=12)
    for u,v,data in G.edges(data=True):
        w=data["weight"]; x=[pos[u][0],pos[v][0]]; y=[pos[u][1],pos[v][1]]
        ax.plot(x,y,color="#aaaaaa",alpha=min(0.08+w*0.07,0.6),lw=min(w*0.6,3),zorder=1)
    for acc in accs:
        dom=dom_bp(acc); col=tcmap.get(dom,"#888888"); x,y=pos[acc]
        ax.scatter(x,y,s=300+G.degree(acc)*80,color=col,alpha=0.92,edgecolors="white",linewidths=1.2,zorder=3)
        ax.text(x,y+0.045,G.nodes[acc]["label"],ha="center",va="bottom",fontsize=7,color="white",
                fontweight="bold",zorder=4,bbox=dict(boxstyle="round,pad=0.15",facecolor="#111",edgecolor="none",alpha=0.65))
    ax.legend(handles=[mpatches.Patch(color=tcmap[t],label=t) for t in top_t]+
              [mpatches.Patch(color="#888888",label="other")],
              loc="lower right",fontsize=9,framealpha=0.8,facecolor="#111",labelcolor="white",
              title="Dominant GO-BP",title_fontsize=9)
    ax.axis("off"); plt.tight_layout(); savefig("Fig29_PPI_Network_CoreProteome.png",fig)
    node_rows=[{"accession":a,"label":G.nodes[a]["label"],"degree":G.degree(a),
                "dominant_GO_BP":dom_bp(a)} for a in accs]
    edge_rows=[{"protein_a":u,"protein_b":v,"weight":d["weight"]} for u,v,d in G.edges(data=True)]
    save_csv("Fig29_PPI_Network_CoreProteome_nodes.csv", pd.DataFrame(node_rows))
    save_csv("Fig29_PPI_Network_CoreProteome_edges.csv", pd.DataFrame(edge_rows))

def fig30_chord_pathways(enrich_results, keys, shorts, colors):
    pw_hits={}; all_pw=set()
    for ki,k in enumerate(keys):
        for db in ["REACTOME","KEGG","WIKIPATHWAYS"]:
            res=enrich_results[k][db]
            for _,row in res.iterrows():
                ps=pw_short(row["pathway"]); all_pw.add(ps)
                pw_hits[(ki,ps)]=max(pw_hits.get((ki,ps),0),int(row["study_hits"]))
    def pw_score(p): return (-sum(1 for ki in range(4) if pw_hits.get((ki,p),0)>0),
                             -sum(pw_hits.get((ki,p),0) for ki in range(4)))
    top_pw=sorted(all_pw,key=pw_score)[:16]
    if not top_pw: print("  [Fig30] No pathways."); return
    n_sp=4; n_pw=len(top_pw); n_nodes=n_sp+n_pw
    flow_mat=np.zeros((n_nodes,n_nodes))
    for ki in range(4):
        for pi,pw in enumerate(top_pw):
            w=pw_hits.get((ki,pw),0)
            if w>0: flow_mat[ki][n_sp+pi]=w; flow_mat[n_sp+pi][ki]=w
    node_totals=flow_mat.sum(axis=1)
    if node_totals.sum()==0: print("  [Fig30] No flow data."); return
    gap_deg=3; total_arc=360-n_nodes*gap_deg
    angles={}; cum=0
    for i in range(n_nodes):
        ext=max(node_totals[i]/node_totals.sum()*total_arc,2)
        angles[i]=(cum+i*gap_deg,cum+i*gap_deg+ext); cum+=ext
    node_colors=list(colors)+[tuple(plt.cm.Set3(j/n_pw)) for j in range(n_pw)]
    def d2r(d): return d*math.pi/180
    fig,ax=plt.subplots(figsize=(16,16),subplot_kw={"aspect":"equal"})
    ax.set_facecolor("#0a0a1a"); fig.patch.set_facecolor("#0a0a1a")
    ax.set_xlim(-2.3,2.3); ax.set_ylim(-2.3,2.3); ax.axis("off")
    ax.set_title("Figure 30: Comparative Chord Diagram – Coral Species × Enriched Pathways\n"
                 "(Arc width = protein count; ribbons = enrichment associations)",
                 fontsize=12,fontweight="bold",color="white",pad=14)
    arc_r=1.0
    for i in range(n_nodes):
        a1,a2=angles[i]; thetas=np.linspace(d2r(a1),d2r(a2),100)
        ax.plot(arc_r*np.cos(thetas),arc_r*np.sin(thetas),
                color=node_colors[i],lw=14,alpha=0.92,solid_capstyle="round",zorder=5)
        mid=(a1+a2)/2; lbl=(f"$\\it{{{shorts[i]}}}$" if i<4 else top_pw[i-4][:32])
        ax.text(1.22*math.cos(d2r(mid)),1.22*math.sin(d2r(mid)),lbl,
                ha="center",va="center",fontsize=9 if i<4 else 7.5,color=node_colors[i],
                fontweight="bold" if i<4 else "normal",
                bbox=dict(boxstyle="round,pad=0.2",facecolor="#0a0a1a",
                          edgecolor=node_colors[i] if i<4 else "none",alpha=0.85))
    offsets={i:angles[i][0] for i in range(n_nodes)}
    extents={i:(angles[i][1]-angles[i][0]) for i in range(n_nodes)}
    for ki in range(4):
        for pi in range(n_pw):
            ni=n_sp+pi; w=flow_mat[ki][ni]
            if w==0 or node_totals[ki]==0 or node_totals[ni]==0: continue
            ei=w/node_totals[ki]*extents[ki]; ej=w/node_totals[ni]*extents[ni]
            ais=offsets[ki]; aie=ais+ei; ajs=offsets[ni]; aje=ajs+ej
            offsets[ki]+=ei; offsets[ni]+=ej
            t=np.linspace(0,1,80); O=np.array([0.0,0.0])
            p0=np.array([arc_r*math.cos(d2r(ais)),arc_r*math.sin(d2r(ais))])
            p1=np.array([arc_r*math.cos(d2r(aie)),arc_r*math.sin(d2r(aie))])
            p2=np.array([arc_r*math.cos(d2r(ajs)),arc_r*math.sin(d2r(ajs))])
            p3=np.array([arc_r*math.cos(d2r(aje)),arc_r*math.sin(d2r(aje))])
            xs1=(1-t)**2*p0[0]+2*(1-t)*t*O[0]+t**2*p3[0]; ys1=(1-t)**2*p0[1]+2*(1-t)*t*O[1]+t**2*p3[1]
            xs2=(1-t)**2*p1[0]+2*(1-t)*t*O[0]+t**2*p2[0]; ys2=(1-t)**2*p1[1]+2*(1-t)*t*O[1]+t**2*p2[1]
            ax.fill(np.concatenate([xs1,xs2[::-1]]),np.concatenate([ys1,ys2[::-1]]),
                    color=colors[ki],alpha=0.30,zorder=1)
    ax.text(0,0,"Species\n×\nPathways",ha="center",va="center",fontsize=10,color="white",fontweight="bold",alpha=0.5)
    ax.legend(handles=[mpatches.Patch(color=colors[i],label=f"$\\it{{{shorts[i]}}}$") for i in range(4)],
              loc="lower right",fontsize=9,facecolor="#111",labelcolor="white",framealpha=0.8,title="Species",title_fontsize=9)
    plt.tight_layout(); savefig("Fig30_Chord_Species_Pathways.png",fig)
    csv_rows=[{"species":shorts[ki],"pathway":pw,"protein_hits":pw_hits.get((ki,pw),0)}
              for ki in range(4) for pw in top_pw if pw_hits.get((ki,pw),0)>0]
    save_csv("Fig30_Chord_Species_Pathways.csv", pd.DataFrame(csv_rows))

def fig31_ppi_all(dfs, keys, shorts, colors, pw_mem):
    all_proteins=pd.concat(list(dfs.values()),ignore_index=True).drop_duplicates(subset=COL_ACCESSION)
    sets_acc=[set(dfs[k][COL_ACCESSION]) for k in keys]
    def presence(acc): return sum(1 for s in sets_acc if acc in s)
    all_proteins["presence"]=all_proteins[COL_ACCESSION].apply(presence)
    all_proteins["score"]=all_proteins[COL_PEPTIDES]*all_proteins[COL_COVERAGE]*all_proteins["presence"]
    top40=all_proteins.nlargest(40,"score").copy()
    top40_names={r[COL_ACCESSION]:pname(r[COL_DESC]) for _,r in top40.iterrows()}
    top40_pws={acc:pw_mem.get(acc, set()) for acc in top40_names}
    top40_gos={r[COL_ACCESSION]:get_go_set(r) for _,r in top40.iterrows()}
    G=nx.Graph()
    for _,row in top40.iterrows():
        acc=row[COL_ACCESSION]
        G.add_node(acc,label=top40_names[acc][:22],pw=top40_pws[acc],
                   go=top40_gos.get(acc,set()),presence=int(row["presence"]),peptides=int(row[COL_PEPTIDES]))
    accs=list(top40[COL_ACCESSION])
    for i in range(len(accs)):
        for j in range(i+1,len(accs)):
            a,b=accs[i],accs[j]
            w=len(G.nodes[a]["pw"]&G.nodes[b]["pw"])*2+len(G.nodes[a]["go"]&G.nodes[b]["go"])
            if w>0: G.add_edge(a,b,weight=w)
    print(f"  PPI all: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    presence_color={1:colors[0],2:colors[1],3:colors[2],4:"#C0392B"}
    pos=nx.spring_layout(G,k=2.4,seed=7,weight="weight",iterations=100)
    fig,ax=plt.subplots(figsize=(20,18)); ax.set_facecolor("#0d0d1a"); fig.patch.set_facecolor("#0d0d1a")
    ax.set_title("Figure 31: STRING-style PPI Co-association Network – Top 40 Biologically Important Proteins\n"
                 "(Node colour = # species; node size = peptide count; edge width = shared annotations)",
                 fontsize=12,fontweight="bold",color="white",pad=14)
    max_w=max((d["weight"] for _,_,d in G.edges(data=True)),default=1)
    for u,v,data in G.edges(data=True):
        w=data["weight"]; x=[pos[u][0],pos[v][0]]; y=[pos[u][1],pos[v][1]]
        ax.plot(x,y,color="#aaaaaa",alpha=min(0.08+w*0.06,0.5),lw=min(w*0.5,3),zorder=1)
    for acc in accs:
        node=G.nodes[acc]; pres=node["presence"]; pep=node["peptides"]
        col=presence_color.get(pres,"#888"); x,y=pos[acc]
        ax.scatter(x,y,s=100+pep*20,color=col,alpha=0.88,edgecolors="white",linewidths=1.2,zorder=3)
        ax.text(x,y+0.04,node["label"],ha="center",va="bottom",fontsize=6.5,color="white",
                fontweight="bold",zorder=4,bbox=dict(boxstyle="round,pad=0.12",facecolor="#111",edgecolor="none",alpha=0.6))
    for pres,lbl in [(1,"1 species"),(2,"2 species"),(3,"3 species"),(4,"4 species (core)")]:
        ax.scatter([],[],s=200,color=presence_color[pres],label=lbl,edgecolors="white",linewidths=1)
    ax.legend(fontsize=9,loc="lower left",facecolor="#111",labelcolor="white",
              framealpha=0.8,title="Node: # species",title_fontsize=9)
    ax.axis("off"); plt.tight_layout(); savefig("Fig31_PPI_Network_AllSpecies.png",fig)
    node_rows=[{"accession":acc,"label":G.nodes[acc]["label"],"n_species_present":G.nodes[acc]["presence"],
                "peptides":G.nodes[acc]["peptides"],"degree":G.degree(acc)} for acc in accs]
    edge_rows=[{"protein_a":u,"protein_b":v,"weight":d["weight"]} for u,v,d in G.edges(data=True)]
    save_csv("Fig31_PPI_Network_AllSpecies_nodes.csv", pd.DataFrame(node_rows))
    save_csv("Fig31_PPI_Network_AllSpecies_edges.csv", pd.DataFrame(edge_rows))

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*70)
    print("  CORAL PROTEOMICS FULL PIPELINE  (Figures 1–31)")
    print("  Gulf of Mannar | BRIC-RGCB × VIT Vellore")
    print("="*70)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    print(f"\n[1/5] Output folders: ./{OUT_DIR}/  (figures)  +  ./{CSV_DIR}/  (csv data)\n")

    print("[2/5] Loading data …")
    dfs = load_samples()
    keys   = list(SAMPLE_FILES.keys())
    colors = [SAMPLE_FILES[k]["color"] for k in keys]
    shorts = [sn(SAMPLE_FILES[k]["species"]) for k in keys]
    print()

    print("[3/5] Fetching real GO + pathway annotations (UniProt/KEGG/Reactome/WikiPathways) …")
    dfs, reactome_mem, kegg_mem, wp_mem, coverage = build_real_annotations(dfs)
    save_csv("Annotation_Fetch_Coverage.csv", pd.DataFrame(coverage))
    top_go_terms_summary(dfs, keys, shorts)
    pw_mem = {}  # merged Reactome+KEGG+WikiPathways membership, used for PPI network edge-weighting
    for a in coverage["accession"]:
        pw_mem[a] = reactome_mem.get(a, set()) | kegg_mem.get(a, set()) | wp_mem.get(a, set())
    print()

    print("[4/5] Running enrichment analyses …")
    enrich, uniq_enrich, core_enrich, unique_dfs, core_df, bg = run_all_enrichment(
        dfs, keys, reactome_mem, kegg_mem, wp_mem)
    print()

    print("[5/5] Generating 31 figures …")
    # Descriptive / comparative (1–19)
    fig1(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig2(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig3(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig4(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig5(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig6(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig7(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig8(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig9(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig10(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig11(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig12(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig13(dfs, SAMPLE_FILES, keys, shorts, colors, reactome_mem, kegg_mem, wp_mem)
    fig14(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig15(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig16(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig17(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig18(dfs, SAMPLE_FILES, keys, shorts, colors)
    fig19(SAMPLE_FILES, keys, shorts, colors)
    # Enrichment (20–27)
    fig_enrichment_bars(enrich, keys, shorts, colors, "REACTOME",     20, "Reactome",     dfs)
    fig_enrichment_bars(enrich, keys, shorts, colors, "WIKIPATHWAYS", 21, "WikiPathways", dfs)
    fig_enrichment_bars(enrich, keys, shorts, colors, "KEGG",         22, "KEGG",         dfs)
    fig_comparative_heatmap(enrich, keys, shorts, colors, "REACTOME",     23, "Reactome")
    fig_comparative_heatmap(enrich, keys, shorts, colors, "KEGG",         24, "KEGG")
    fig_comparative_heatmap(enrich, keys, shorts, colors, "WIKIPATHWAYS", 25, "WikiPathways")
    fig26_unique_go(dfs, unique_dfs, keys, shorts, colors)
    fig27_unique_pathway(uniq_enrich, unique_dfs, keys, shorts, colors)
    # Core + PPI + Chord (28–31)
    fig28_core_enrichment(core_enrich, core_df)
    fig29_ppi_core(core_df, pw_mem)
    fig30_chord_pathways(enrich, keys, shorts, colors)
    fig31_ppi_all(dfs, keys, shorts, colors, pw_mem)

    print(f"\n{'='*70}")
    print(f"  Done! 31 figures saved to → ./{OUT_DIR}/")
    print(f"  Matching CSV data saved to → ./{CSV_DIR}/")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
