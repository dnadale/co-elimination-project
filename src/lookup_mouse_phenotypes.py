"""Look up mouse orthologs and knockout phenotypes (IMPC + MGI via Alliance) for the 19 unlinked-pair genes."""
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "vetting"
OUT.mkdir(parents=True, exist_ok=True)


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)
    except Exception:
        return None


pairs = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t")
sym = pd.read_csv(ROOT / "raw_data/ensembl/ensembl_id_symbol_map.tsv", sep="\t").drop_duplicates("Gene stable ID").set_index("Gene stable ID")
ids = list(dict.fromkeys(list(pairs.gene_a) + list(pairs.gene_b)))

rows = []
for gid in ids:
    g0 = gid.split(".")[0]
    name = sym["Gene name"].get(g0)
    desc = str(sym["Gene description"].get(g0, "")).split(" [")[0]
    mouse_sym = mouse_id = mgi = None
    d = get(f"https://rest.ensembl.org/homology/id/human/{g0}?target_species=mus_musculus;type=orthologues;content-type=application/json")
    homs = (((d or {}).get("data") or [{}])[0]).get("homologies", [])
    if homs:
        mouse_id = homs[0]["target"]["id"]
        x = get(f"https://rest.ensembl.org/xrefs/id/{mouse_id}?external_db=MGI;content-type=application/json") or []
        if x:
            mgi = x[0]["primary_id"]
            mouse_sym = x[0]["display_id"]
    impc, mgi_terms = [], []
    if mouse_sym:
        q = urllib.parse.quote(f"marker_symbol:{mouse_sym}")
        s = get(f"https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q={q}&rows=200&wt=json")
        docs = ((s or {}).get("response") or {}).get("docs", [])
        impc = sorted({f"{d['mp_term_name']} ({d.get('zygosity', '?')})" for d in docs})
    if mgi:
        a = get(f"https://www.alliancegenome.org/api/gene/{mgi}/phenotypes?limit=200")
        for r in (a or {}).get("results", []):
            if r.get("phenotypeStatement"):
                mgi_terms.append(r["phenotypeStatement"])
        mgi_terms = sorted(set(mgi_terms))
    rows.append(dict(
        human_gene=name or g0, ensembl_id=g0, description=desc,
        mouse_ortholog=mouse_sym or "NONE", mgi_id=mgi or "",
        n_impc_phenotypes=len(impc), impc_phenotypes="; ".join(impc[:15]),
        n_mgi_phenotypes=len(mgi_terms), mgi_phenotypes="; ".join(mgi_terms[:15]),
    ))
    print(rows[-1]["human_gene"], rows[-1]["mouse_ortholog"], len(impc), len(mgi_terms), flush=True)

pd.DataFrame(rows).to_csv(OUT / "mouse_knockout_phenotypes.tsv", sep="\t", index=False)
