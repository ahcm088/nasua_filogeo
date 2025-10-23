#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fst_nasua.py
Calcula Hudson's F_ST par-a-par (mtDNA, haploide) para Nasua nasua e Nasua narica,
separadamente, a partir de um FASTA ALINHADO e um CSV de metadados.

Estrutura esperada:
- dados/
    - nasua_alinhado.fasta
    - metadata.csv   (cabeçalho: Accession,Organism,Country,Geo_loc_name,Collector,Collection_date,Locality,Sequence_length,Definition)
- resultados/        (será criado)

Uso:
    python scripts/fst_nasua.py \
        --fasta dados/nasua_alinhado.fasta \
        --metadata dados/metadata.csv \
        --pop-field Geo_loc_name \
        --min-n 3 \
        --resultados resultados

Autor: você :)
"""

import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import allel
import matplotlib.pyplot as plt

SPECIES_LIST = ["Nasua nasua", "Nasua narica"]  # ajuste se precisar

# ---------------------- utilidades ----------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def read_fasta_alignment(fasta_path):
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    if not recs:
        raise SystemExit(f"[erro] FASTA vazio: {fasta_path}")
    L0 = len(recs[0].seq)
    for r in recs:
        if len(r.seq) != L0:
            raise SystemExit(f"[erro] Alinhamento inconsistente: {r.id} tem {len(r.seq)} bp; esperado {L0} bp.")
    return recs, L0

def seqs_to_array(records):
    # n_seqs x n_sites, caracteres em maiúsculas
    return np.array([list(str(r.seq).upper()) for r in records], dtype='<U1')

def accession_from_header(h):
    # "MK144305.1_Nasua_nasua" -> "MK144305.1"
    return h.split('_')[0]

def species_from_header(h):
    # "MK144305.1_Nasua_nasua" -> "Nasua nasua"
    if '_' in h:
        sp = h.split('_', 1)[1].replace('_', ' ')
        return sp
    return ""

def filter_biallelic_snps(char_array):
    """Mantém apenas sites bialélicos A/C/G/T, ignora Ns/gaps."""
    valid_nt = np.isin(char_array, np.array(list("ACGT")))
    # Sites com pelo menos 2 observações válidas
    good_site = valid_nt.sum(axis=0) >= 2
    arr2 = char_array[:, good_site]
    valid2 = valid_nt[:, good_site]

    keep_cols = []
    for j in range(arr2.shape[1]):
        alleles = np.unique(arr2[valid2[:, j], j])
        if alleles.size == 2:
            keep_cols.append(j)
    keep_cols = np.array(keep_cols, dtype=int)
    return arr2[:, keep_cols], valid2[:, keep_cols]

def to_012_haploid(arr_bi, valid_bi):
    """
    Transforma em genótipo haploide 0/1 com -1 = missing.
    Define ref/alt alfabeticamente por coluna para consistência.
    """
    n, m = arr_bi.shape
    gt = np.full((n, m), -1, dtype=np.int8)
    ref_alt = []
    for j in range(m):
        alleles = np.unique(arr_bi[valid_bi[:, j], j])
        ref, alt = alleles[0], alleles[1]
        ref_alt.append((ref, alt))
        col = arr_bi[:, j]
        gt[(col == ref) & valid_bi[:, j], j] = 0
        gt[(col == alt) & valid_bi[:, j], j] = 1
    return gt, ref_alt

def hudson_pairwise_fst(gt_01, idx_a, idx_b):
    """
    Hudson's F_ST (scikit-allel) entre duas populações.
    gt_01: (n_amostras x n_sites) em {0,1,-1}, haploide
    """
    m = gt_01.shape[1]
    ac1 = np.zeros((2, m), dtype=int)
    ac2 = np.zeros((2, m), dtype=int)

    g1 = gt_01[idx_a, :]
    g2 = gt_01[idx_b, :]

    for j in range(m):
        col1 = g1[:, j]
        ac1[0, j] = np.sum(col1 == 0)
        ac1[1, j] = np.sum(col1 == 1)

        col2 = g2[:, j]
        ac2[0, j] = np.sum(col2 == 0)
        ac2[1, j] = np.sum(col2 == 1)

    num, den = allel.hudson_fst(ac1, ac2)  # por-site
    mask = den > 0
    fst = np.sum(num[mask]) / np.sum(den[mask]) if np.any(mask) else np.nan
    return fst, int(np.sum(num[mask])), int(np.sum(den[mask])), int(mask.sum())

def plot_heatmap(matrix_df, out_png, title="Hudson F_ST (pairwise)"):
    plt.figure(figsize=(7, 6))
    # Ordenar e desenhar
    pops = list(matrix_df.index)
    M = matrix_df.values.astype(float)
    im = plt.imshow(M, interpolation='nearest')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(range(len(pops)), pops, rotation=45, ha='right')
    plt.yticks(range(len(pops)), pops)
    plt.title(title)
    # Anotar valores
    for i in range(len(pops)):
        for j in range(len(pops)):
            val = M[i, j]
            txt = "" if (i == j or not np.isfinite(val)) else f"{val:.3f}"
            plt.text(j, i, txt, ha='center', va='center', fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# ---------------------- pipeline ----------------------
def run_species(species_name, records, meta_df, pop_field, out_dir, min_n=3):
    print(f"\n===> Espécie: {species_name}")

    # Tabela auxiliar a partir do FASTA
    headers = [r.id for r in records]
    df_fa = pd.DataFrame({
        "header": headers,
        "Accession": [accession_from_header(h) for h in headers],
        "Organism_from_header": [species_from_header(h) for h in headers],
    })

    # Merge com metadata por Accession
    df = df_fa.merge(meta_df, on="Accession", how="left")

    # Filtra espécie pelo campo Organism do CSV
    df_sp = df[df["Organism"].astype(str).str.strip() == species_name].copy()
    if df_sp.empty:
        print(f"[aviso] Sem sequências para {species_name} após o merge por Accession.")
        return

    # Reordenar os registros na ordem de df_sp
    id_to_rec = {r.id: r for r in records}
    rec_sp = [id_to_rec[h] for h in df_sp["header"]]
    arr = seqs_to_array(rec_sp)

    # Descobrir SNPs bialélicos
    arr_bi, valid_bi = filter_biallelic_snps(arr)
    if arr_bi.size == 0:
        print("[erro] Nenhum SNP bialélico após filtragem (A/C/G/T). Verifique alinhamento e qualidade.")
        return
    gt, _ = to_012_haploid(arr_bi, valid_bi)

    # Definir populações
    if pop_field not in df_sp.columns:
        raise SystemExit(f"[erro] Campo '{pop_field}' não existe no metadata.csv.")
    pops_series = df_sp[pop_field].astype(str).str.strip().replace({"": np.nan, "nan": np.nan})
    pop_counts = pops_series.value_counts(dropna=True)
    keep_pops = pop_counts[pop_counts >= min_n].index.tolist()
    if len(keep_pops) < 2:
        print(f"[aviso] Menos de duas populações com n >= {min_n}. Mantidas: {keep_pops}")
        return
    pops_idx = {p: np.where(pops_series.values == p)[0] for p in keep_pops}

    # Salvar tamanhos amostrais
    pop_sizes = pd.DataFrame({"populacao": keep_pops, "n": [len(pops_idx[p]) for p in keep_pops]})
    pop_sizes.to_csv(os.path.join(out_dir, "pop_sizes.csv"), index=False)

    # Matriz F_ST par-a-par
    pops = keep_pops
    nP = len(pops)
    M = pd.DataFrame(np.nan, index=pops, columns=pops, dtype=float)
    details = []
    for i in range(nP):
        for j in range(i + 1, nP):
            p1, p2 = pops[i], pops[j]
            fst, num, den, n_sites = hudson_pairwise_fst(gt, pops_idx[p1], pops_idx[p2])
            M.loc[p1, p2] = fst
            M.loc[p2, p1] = fst
            details.append({
                "pop1": p1, "pop2": p2,
                "fst_hudson": fst,
                "num_sum": num, "den_sum": den,
                "n_sites": n_sites,
                "n_pop1": len(pops_idx[p1]), "n_pop2": len(pops_idx[p2])
            })

    # Salvar resultados
    M.to_csv(os.path.join(out_dir, "fst_pairwise_hudson.csv"))
    pd.DataFrame(details).to_csv(os.path.join(out_dir, "fst_pair_details.csv"), index=False)
    print(f"[ok] F_ST salvo em: {os.path.join(out_dir,'fst_pairwise_hudson.csv')}")
    print(f"[ok] Detalhes em:  {os.path.join(out_dir,'fst_pair_details.csv')}")

    # Heatmap
    heat_png = os.path.join(out_dir, "fst_pairwise_hudson_heatmap.png")
    plot_heatmap(M, heat_png, title=f"Hudson F_ST — {species_name}")
    print(f"[ok] Heatmap salvo em: {heat_png}")

def main():
    ap = argparse.ArgumentParser(description="Hudson's F_ST pairwise para Nasua nasua e Nasua narica (mtDNA haploide).")
    ap.add_argument("--fasta", default="dados/nasua_alinhado.fasta", help="FASTA ALINHADO com cabeçalhos tipo >Accession_Especie")
    ap.add_argument("--metadata", default="dados/metadata.csv", help="CSV de metadados")
    ap.add_argument("--resultados", default="resultados", help="Pasta de saída")
    ap.add_argument("--pop-field", default="Geo_loc_name", help="Campo do CSV que define populações (ex.: Geo_loc_name, Country, Locality)")
    ap.add_argument("--min-n", type=int, default=3, help="Amostras mínimas por população")
    args = ap.parse_args()

    # Checagens e leitura
    if not os.path.exists(args.fasta):
        raise SystemExit(f"[erro] FASTA não encontrado: {args.fasta}")
    if not os.path.exists(args.metadata):
        raise SystemExit(f"[erro] CSV não encontrado: {args.metadata}")
    ensure_dir(args.resultados)

    records, L = read_fasta_alignment(args.fasta)
    print(f"[info] FASTA: {args.fasta} | {len(records)} seqs | {L} bp alinhados")
    meta = pd.read_csv(args.metadata)
    print(f"[info] CSV:   {args.metadata} | {len(meta)} linhas")

    # Rodar espécies
    for sp in SPECIES_LIST:
        out_dir_sp = os.path.join(args.resultados, sp.replace(" ", "_"))
        ensure_dir(out_dir_sp)
        run_species(sp, records, meta, args.pop_field, out_dir_sp, args.min_n)

    print("\n[feito] Consulte a pasta de resultados.")

if __name__ == "__main__":
    main()
