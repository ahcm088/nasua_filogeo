#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, argparse, numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from itertools import combinations

def acc_from_header(h): return h.split('_')[0]
def sp_from_header(h):  return h.split('_',1)[1].replace('_',' ') if '_' in h else ""

def read_fasta_aligned(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs: raise SystemExit(f"[erro] FASTA vazio: {path}")
    L = len(recs[0].seq)
    for r in recs:
        if len(r.seq) != L:
            raise SystemExit("[erro] FASTA não alinhado (comprimentos distintos).")
    return recs, L

def pairwise_diff(seqA, seqB):
    dif=0; L=0
    for a,b in zip(seqA, seqB):
        if a in "ACGT" and b in "ACGT":
            L += 1
            if a!=b: dif += 1
    return dif, L

def avg_pairwise_within(seqs):
    # retorna K (diferença média) e Lmédio (bases comparadas)
    if len(seqs) < 2: return (np.nan, 0)
    tot_d=0; tot_L=0; m=0
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            d,L = pairwise_diff(seqs[i], seqs[j])
            if L>0:
                tot_d += d; tot_L += L; m += 1
    if m==0: return (np.nan, 0)
    return (tot_d/m, tot_L/m)

def tajimas_D(n, S, pi, L):
    """
    n: número de amostras na população
    S: nº de sítios segregantes
    pi: nucleotide diversity (diferença média por sítio)
    L: nº de sítios comparáveis (use L_effective se calculado)
    Retorna (D, D) ou (NaN, NaN) quando indefinido.
    """
    # condições que inviabilizam o cálculo
    if n is None or n < 2 or L is None or L <= 0 or not np.isfinite(pi) or S is None or S <= 0:
        return np.nan, np.nan

    # coeficientes de Tajima (1989)
    a1 = np.sum(1.0 / np.arange(1, n))
    a2 = np.sum(1.0 / (np.arange(1, n) ** 2))
    if a1 == 0:
        return np.nan, np.nan

    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n**2 + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    thetaW = (S / a1) / L  # por sítio
    var_term = e1 * S + e2 * S * (S - 1.0)

    # guarda contra zeros/negativos e precisão numérica
    if not np.isfinite(thetaW) or var_term <= 0:
        return np.nan, np.nan

    with np.errstate(invalid="ignore", divide="ignore"):
        D = (pi - thetaW) / np.sqrt(var_term)
    return (D, D)

def haplotype_diversity(hap_ids):
    """
    Hd = n/(n-1) * (1 - sum(p_i^2)), Nei (1987)
    """
    n = len(hap_ids)
    if n <= 1: return np.nan
    vals, counts = np.unique(hap_ids, return_counts=True)
    p2 = np.sum((counts/n)**2)
    return (n/(n-1.0))*(1.0 - p2)

def run_for_species(species, fasta_path, meta_path, pop_field, outdir):
    os.makedirs(outdir, exist_ok=True)
    recs, L = read_fasta_aligned(fasta_path)
    meta = pd.read_csv(meta_path)

    headers = [r.id for r in recs]
    df = pd.DataFrame({"header": headers,
                       "Accession":[acc_from_header(h) for h in headers],
                       "sp_hdr":[sp_from_header(h) for h in headers]})
    df = df.merge(meta, on="Accession", how="left")
    df_sp = df[df["Organism"].astype(str).str.strip()==species].copy()
    if df_sp.empty:
        print(f"[aviso] Sem sequências para {species}."); return

    id2seq = {r.id:str(r.seq).upper() for r in recs}
    df_sp["seq"] = df_sp["header"].map(id2seq)
    df_sp[pop_field] = df_sp[pop_field].astype(str).str.strip()

    # --- métricas por população ---
    rows = []
    for pop, grp in df_sp.groupby(pop_field):
        seqs = grp["seq"].tolist()
        n = len(seqs)
        # K e pi
        K, Lm = avg_pairwise_within(seqs)
        pi = (K / Lm) if (Lm and np.isfinite(K)) else np.nan
        # S (segregating sites) e L_effective
        if n >= 2:
            arr = np.array([list(s) for s in seqs])
            S=0; Luse=0
            for col in arr.T:
                ok = np.isin(col, list("ACGT"))
                alle = np.unique(col[ok])
                if ok.sum() >= 2:
                    Luse += 1
                    if len(alle) >= 2:
                        S += 1
        else:
            S=0; Luse=0
        # haplotype diversity
        hap_ids = pd.Series(seqs).astype('category').cat.codes.astype(str)  # id por sequência idêntica
        Hd = haplotype_diversity(hap_ids)
        # Tajima's D (seguro)
        D, z = tajimas_D(n=n, S=S, pi=pi, L=Luse if Luse>0 else L)

        rows.append({"species": species, "population": pop, "n": n,
                     "n_hap": int(pd.Series(hap_ids).nunique()),
                     "Hd": Hd, "pi": pi, "K": K, "S": S, "L_effective": Luse if Luse>0 else L,
                     "TajimaD": D})
    bypop = pd.DataFrame(rows).sort_values(["species","population"])
    bypop.to_csv(os.path.join(outdir, "diversity_by_population.csv"), index=False)

    # --- Dxy/K entre pares de populações ---
    pairs = []
    pops = bypop["population"].tolist()
    seqs_by_pop = {p: df_sp[df_sp[pop_field]==p]["seq"].tolist() for p in pops}
    for p1, p2 in combinations(pops, 2):
        A = seqs_by_pop[p1]; B = seqs_by_pop[p2]
        if len(A)==0 or len(B)==0: continue
        tot_d=0; tot_L=0; m=0
        for a in A:
            for b in B:
                d,L = pairwise_diff(a,b)
                if L>0:
                    tot_d += d; tot_L += L; m += 1
        if m==0:
            Kxy=np.nan; Dxy=np.nan
        else:
            Kxy = tot_d/m
            Dxy = Kxy / (tot_L/m)
        pairs.append({"species":species, "pop1":p1, "pop2":p2, "Kxy":Kxy, "Dxy":Dxy})
    bypair = pd.DataFrame(pairs)
    bypair.to_csv(os.path.join(outdir, "dxy_pairwise.csv"), index=False)

    # --- plots simples: pi e Hd ---
    def barplot(col, fname, ylabel):
        dfp = bypop.dropna(subset=[col])
        if dfp.empty: return
        plt.figure(figsize=(7,4))
        plt.bar(dfp["population"], dfp[col])
        plt.xticks(rotation=45, ha='right'); plt.ylabel(ylabel)
        plt.title(f"{species} — {ylabel} por população")
        plt.tight_layout(); plt.savefig(os.path.join(outdir,fname), dpi=300); plt.close()
    barplot("pi", "pi_by_population.png", "π (por sítio)")
    barplot("Hd", "Hd_by_population.png", "Haplotype diversity (Hd)")

    print(f"[ok] {species}: métricas salvas em {outdir}")

def main():
    ap = argparse.ArgumentParser(description="Diversidade por população (π, Hd, S, TajimaD) e Dxy entre populações.")
    ap.add_argument("--fasta", default="dados/Nasua_F_ST_aligned.fasta")
    ap.add_argument("--metadata", default="dados/Nasua_CYTB_metadata_com_localizacao.csv")
    ap.add_argument("--pop-field", default="Geo_loc_name")
    ap.add_argument("--resultados", default="resultados")
    ap.add_argument("--species", choices=["Nasua nasua","Nasua narica","both"], default="both")
    args = ap.parse_args()

    if args.species in ("Nasua nasua","Nasua narica"):
        sp_list=[args.species]
    else:
        sp_list=["Nasua nasua","Nasua narica"]

    for sp in sp_list:
        outdir = os.path.join(args.resultados, sp.replace(" ","_"))
        run_for_species(sp, args.fasta, args.metadata, args.pop_field, outdir)

if __name__ == "__main__":
    main()
