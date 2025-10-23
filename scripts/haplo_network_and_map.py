#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, argparse, numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter, defaultdict

def acc_from_header(h):  # "MK144305.1_Nasua_nasua" -> "MK144305.1"
    return h.split('_')[0]

def species_from_header(h):  # -> "Nasua nasua"
    return h.split('_', 1)[1].replace('_', ' ') if '_' in h else ""

def read_fasta_aligned(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        raise SystemExit(f"[erro] FASTA vazio: {path}")
    L = len(recs[0].seq)
    for r in recs:
        if len(r.seq) != L:
            raise SystemExit("FASTA não está alinhado (comprimentos diferentes).")
    return recs, L

def hamming(a, b):
    return sum(x!=y for x,y in zip(a,b))

def build_haplotypes(records, df_meta, species, pop_field):
    # filtra espécie via metadata (coluna "Organism")
    headers = [r.id for r in records]
    df_fa = pd.DataFrame({"header": headers,
                          "Accession":[acc_from_header(h) for h in headers],
                          "Org_hdr":[species_from_header(h) for h in headers]})
    df = df_fa.merge(df_meta, on="Accession", how="left")
    df_sp = df[df["Organism"].astype(str).str.strip()==species].copy()
    if df_sp.empty:
        raise SystemExit(f"[erro] Sem sequências para {species} após merge por Accession.")
    # monta lista de sequências na ordem de df_sp
    id2rec = {r.id:r for r in records}
    recs_sp = [id2rec[h] for h in df_sp["header"]]
    seqs = [str(r.seq).upper() for r in recs_sp]

    # agrupar por sequência (haplótipos)
    hap_dict = {}   # seq -> hap_id
    hap_list = []   # lista de sequências únicas
    for s in seqs:
        if s not in hap_dict:
            hap_dict[s] = f"H{len(hap_list)+1}"
            hap_list.append(s)
    # mapeamento indivíduo -> hap
    assignment = []
    for h, row in df_sp.iterrows():
        seq = str(id2rec[row["header"]].seq).upper()
        assignment.append({
            "Accession": row["Accession"],
            "Organism": row["Organism"],
            "Population": str(row[pop_field]),
            "Haplotype": hap_dict[seq]
        })
    assign_df = pd.DataFrame(assignment)

    # tamanho de cada hap (n amostras)
    counts = assign_df["Haplotype"].value_counts().to_dict()

    # composição por população (para colorir nós)
    comp = assign_df.groupby(["Haplotype","Population"]).size().unstack(fill_value=0)

    # matriz de distâncias entre haplótipos (Hamming em alinhado)
    H = len(hap_list)
    dist = np.zeros((H,H), dtype=int)
    for i in range(H):
        for j in range(i+1,H):
            d = hamming(hap_list[i], hap_list[j])
            dist[i,j]=dist[j,i]=d

    # grafo MST
    G = nx.Graph()
    hap_ids = [f"H{i+1}" for i in range(H)]
    for i,hid in enumerate(hap_ids):
        G.add_node(hid, size=counts.get(hid,0))
    for i in range(H):
        for j in range(i+1,H):
            if dist[i,j]>0:
                G.add_edge(hap_ids[i], hap_ids[j], weight=dist[i,j])
    if G.number_of_edges()>0:
        T = nx.minimum_spanning_tree(G, weight='weight')
    else:
        T = G.copy()

    # DataFrame com haplótipos
    hap_df = pd.DataFrame({
        "Haplotype": hap_ids,
        "Size": [counts.get(h,0) for h in hap_ids],
        "Sequence": hap_list
    })

    return {
        "hap_df": hap_df,
        "assign_df": assign_df,
        "comp": comp,         # composição por população
        "mst": T,
        "hap_ids": hap_ids
    }

def draw_pie(ax, center, sizes, total, radius=0.25):
    # desenha fatias proporcionais em um círculo (pie "manual")
    if total<=0: return
    import math
    start=0.0
    for frac in sizes:
        if frac<=0: continue
        theta1, theta2 = 2*math.pi*start, 2*math.pi*(start+frac)
        theta = np.linspace(theta1, theta2, 50)
        x = center[0] + radius*np.cos(theta)
        y = center[1] + radius*np.sin(theta)
        ax.fill(np.r_[center[0], x, center[0]], np.r_[center[1], y, center[1]], alpha=0.9)
        start += frac

def plot_haplo_network(res, out_png, title):
    T = res["mst"]
    comp = res["comp"]  # Haplotype x Population
    hap_sizes = res["hap_df"].set_index("Haplotype")["Size"].to_dict()

    # paleta simples por população
    pops = list(comp.columns)
    cmap = plt.cm.get_cmap('tab10', max(10, len(pops)))
    pop2col = {p: cmap(i%10) for i,p in enumerate(pops)}

    pos = nx.spring_layout(T, seed=42)
    plt.figure(figsize=(8,6))

    # desenha arestas com rótulo = peso (diferenças)
    weights = [T[u][v]['weight'] for u,v in T.edges] if T.number_of_edges() else []
    nx.draw_networkx_edges(T, pos, width=[1.5+0.5*(1.0/w) if w>0 else 1 for w in weights])

    # desenha nós como "pizza" (composição por população)
    ax = plt.gca()
    for h in T.nodes:
        # proporções por população
        row = comp.loc[h] if h in comp.index else pd.Series(0,index=pops)
        total = row.sum()
        if total==0:  # sem composição registrada
            ax.scatter([pos[h][0]],[pos[h][1]], s=120, c='grey', zorder=3)
        else:
            fracs = (row/total).values.tolist()
            # set cor ciclicamente para cada fatia
            # truque: mudar a cor ativa do artista antes de cada fill
            start=0
            radius = 0.2 + 0.03*np.log1p(hap_sizes.get(h,1))
            # fatiar manualmente com draw_pie e alternar cor
            import matplotlib as mpl
            colors = [pop2col[p] for p in pops]
            cum=0.0
            for frac, col in zip(fracs, colors):
                if frac<=0: continue
                theta = np.linspace(2*np.pi*cum, 2*np.pi*(cum+frac), 50)
                x = pos[h][0] + radius*np.cos(theta)
                y = pos[h][1] + radius*np.sin(theta)
                ax.fill(np.r_[pos[h][0], x, pos[h][0]], np.r_[pos[h][1], y, pos[h][1]], color=col, alpha=0.95, ec='k', lw=0.4, zorder=3)
                cum += frac
        # borda do nó
        circ = plt.Circle(pos[h], 0.2 + 0.03*np.log1p(hap_sizes.get(h,1)), fill=False, ec='k', lw=0.6, zorder=4)
        ax.add_patch(circ)
        plt.text(pos[h][0], pos[h][1], h, fontsize=8, ha='center', va='center', zorder=5)

    # legenda das populações
    for i,p in enumerate(pops):
        plt.scatter([],[], c=[pop2col[p]], label=p)
    plt.legend(title="População", bbox_to_anchor=(1.04,1), loc="upper left")
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_haplo_map(res, coords_df, out_png, title):
    # composição por população: assign_df -> por localidade, proporção de haplótipos
    assign = res["assign_df"]
    if coords_df is None or coords_df.empty:
        print("[aviso] Sem geo_coords.csv — pulando mapa.")
        return
    # agrega por localidade e haplótipo
    tab = (assign.groupby(["Population","Haplotype"]).size()
           .unstack(fill_value=0))
    df = coords_df.merge(tab, left_on="Geo_loc_name", right_index=True, how="inner")
    if df.empty:
        print("[aviso] Coordenadas não batem com nomes em Population — pulando mapa.")
        return
    plt.figure(figsize=(8,6))
    ax = plt.gca()
    # eixo simples em lon/lat
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    minx, maxx = df["lon"].min()-5, df["lon"].max()+5
    miny, maxy = df["lat"].min()-5, df["lat"].max()+5
    ax.set_xlim(minx, maxx); ax.set_ylim(miny, maxy)
    hap_cols = [c for c in df.columns if c.startswith("H")]
    # paleta por haplótipo
    cmap = plt.cm.get_cmap('tab20', max(20, len(hap_cols)))
    hap2col = {h: cmap(i%20) for i,h in enumerate(hap_cols)}

    for _,row in df.iterrows():
        lat, lon = row["lat"], row["lon"]
        counts = np.array([row[h] for h in hap_cols], dtype=float)
        total = counts.sum()
        if total <= 0:
            continue
        fracs = counts/total
        # desenha pie simples
        start=0.0
        for frac, h in zip(fracs, hap_cols):
            if frac<=0: continue
            theta = np.linspace(2*np.pi*start, 2*np.pi*(start+frac), 50)
            x = lon + 1.0*np.cos(theta)
            y = lat + 1.0*np.sin(theta)
            ax.fill(np.r_[lon, x, lon], np.r_[lat, y, lat], color=hap2col[h], alpha=0.95, ec='k', lw=0.4)
            start += frac
        ax.text(lon, lat, f" {row['Geo_loc_name']} (n={int(total)})", fontsize=8, ha='left', va='center')
    # legenda
    for h in hap_cols:
        plt.scatter([],[], c=[hap2col[h]], label=h)
    plt.legend(title="Haplótipos", bbox_to_anchor=(1.02,1), loc="upper left")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def main():
    ap = argparse.ArgumentParser(description="Rede de haplótipos (MST) e mapa por população.")
    ap.add_argument("--fasta", default="dados/Nasua_F_ST_aligned.fasta")
    ap.add_argument("--metadata", default="dados/Nasua_CYTB_metadata_com_localizacao.csv")
    ap.add_argument("--species", required=True, choices=["Nasua nasua","Nasua narica"])
    ap.add_argument("--pop-field", default="Geo_loc_name")
    ap.add_argument("--coords", default="dados/geo_coords.csv", help="CSV com Geo_loc_name,lat,lon (opcional)")
    ap.add_argument("--outdir", default="resultados")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    species_dir = os.path.join(args.outdir, args.species.replace(" ","_"))
    os.makedirs(species_dir, exist_ok=True)

    recs, L = read_fasta_aligned(args.fasta)
    meta = pd.read_csv(args.metadata)
    coords_df = None
    if args.coords and os.path.exists(args.coords):
        coords_df = pd.read_csv(args.coords)

    res = build_haplotypes(recs, meta, args.species, args.pop_field)

    # salvar CSVs
    res["hap_df"].to_csv(os.path.join(species_dir,"haplotypes.csv"), index=False)
    res["assign_df"].to_csv(os.path.join(species_dir,"haplotype_assignment.csv"), index=False)

    # plots
    plot_haplo_network(res, os.path.join(species_dir,"haplo_network.png"),
                       title=f"Haplótipos CytB — {args.species}")
    plot_haplo_map(res, coords_df, os.path.join(species_dir,"haplo_map.png"),
                   title=f"Haplótipos por localidade — {args.species}")

    print("[ok] Saídas em:", species_dir)

if __name__ == "__main__":
    main()
