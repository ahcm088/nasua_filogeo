#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
from Bio import SeqIO

def accession_from_header(h):
    # Ex.: "MK144305.1_Nasua_nasua" -> "MK144305.1"
    return h.split('_')[0]

def main():
    ap = argparse.ArgumentParser(description="Conta indivíduos por geolocalização (e por espécie x geolocalização) usando FASTA + metadata.")
    ap.add_argument("--fasta", required=True, help="FASTA alinhado (cabecalhos tipo >Accession_Especie)")
    ap.add_argument("--metadata", required=True, help="CSV de metadados (inclui colunas Accession, Organism, Geo_loc_name, ...)")
    ap.add_argument("--outdir", default="resultados", help="Pasta de saída")
    ap.add_argument("--pop-field", default="Geo_loc_name", help="Campo do CSV para populações (ex.: Geo_loc_name, Country, Locality)")
    ap.add_argument("--min-n", type=int, default=3, help="Limite mínimo para destacar populações com n >= min-n")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Acessões presentes no FASTA
    headers = [rec.id for rec in SeqIO.parse(args.fasta, "fasta")]
    if not headers:
        raise SystemExit("[erro] FASTA vazio.")
    accessions = {accession_from_header(h) for h in headers}

    # 2) Metadados
    df = pd.read_csv(args.metadata)
    if "Accession" not in df.columns:
        raise SystemExit("[erro] CSV precisa ter a coluna 'Accession'.")
    if args.pop_field not in df.columns:
        raise SystemExit(f"[erro] Campo '{args.pop_field}' não existe no CSV.")

    # 3) Filtrar somentes os indivíduos cujas acessões estão no FASTA
    df_filt = df[df["Accession"].isin(accessions)].copy()
    if df_filt.empty:
        raise SystemExit("[erro] Após filtrar por Accession do FASTA, não sobraram linhas no CSV.")

    # Limpeza simples do campo de população
    df_filt[args.pop_field] = df_filt[args.pop_field].astype(str).str.strip()
    df_filt.loc[df_filt[args.pop_field].isin(["", "nan", "NA"]), args.pop_field] = pd.NA

    # 4) Contagens
    by_pop = (df_filt
              .groupby(args.pop_field, dropna=True, as_index=False)
              .size()
              .rename(columns={"size": "n"}))
    by_sp_pop = (df_filt
                 .groupby(["Organism", args.pop_field], dropna=True, as_index=False)
                 .size()
                 .rename(columns={"size": "n"}))

    # 5) Salvar
    out1 = os.path.join(args.outdir, "counts_by_geoloc.csv")
    out2 = os.path.join(args.outdir, "counts_by_species_geoloc.csv")
    by_pop.to_csv(out1, index=False)
    by_sp_pop.to_csv(out2, index=False)

    # 6) Mostrar um resumo no console
    print(f"[ok] Salvo: {out1} e {out2}")
    print("\n=== Populações com n >= {0} (todas as espécies juntas) ===".format(args.min_n))
    print(by_pop[by_pop["n"] >= args.min_n].sort_values("n", ascending=False).to_string(index=False))

    print("\n=== Por espécie x população (n >= {0}) ===".format(args.min_n))
    print(by_sp_pop[by_sp_pop["n"] >= args.min_n].sort_values(["Organism","n"], ascending=[True, False]).to_string(index=False))

if __name__ == "__main__":
    main()
