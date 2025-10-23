#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, argparse, numpy as np, pandas as pd
from math import radians, sin, cos, asin, sqrt
import matplotlib.pyplot as plt

# Mantel só será importado quando realmente usado (>=3 pops)
def haversine_km(lat1, lon1, lat2, lon2):
    R=6371.0
    dlat=radians(lat2-lat1); dlon=radians(lon2-lon1)
    a=sin(dlat/2)**2 + cos(radians(lat1))*cos(radians(lat2))*sin(dlon/2)**2
    return 2*R*asin(sqrt(a))

def abbr3(name: str) -> str:
    """Abrevia nome para 3 letras (heurística simples)."""
    if not isinstance(name, str) or not name:
        return name
    # tenta usar iniciais maiúsculas; se não houver, pega 3 primeiras letras
    tokens = [t for t in name.split() if t]
    caps = "".join([t[0] for t in tokens if t[0].isupper()])
    if len(caps) >= 2:
        return (caps[:3]).upper()
    return (name[:3]).upper()

def main():
    ap = argparse.ArgumentParser(description="Mantel F_ST x distância geográfica (com rótulos opcionais e fallback para 2 populações).")
    ap.add_argument("--species-dir", required=True, help="ex.: resultados/Nasua_nasua")
    ap.add_argument("--coords", default="dados/geo_coords.csv", help="CSV: Geo_loc_name,lat,lon")
    ap.add_argument("--label-points", action="store_true", help="Se passado, rotula os pontos do scatter com os pares de países.")
    ap.add_argument("--abbr", action="store_true", help="Se passado, usa abreviações de 3 letras nos rótulos.")
    args = ap.parse_args()

    # --- carregar matriz F_ST ---
    mpath = os.path.join(args.species_dir, "fst_pairwise_hudson.csv")
    if not os.path.exists(mpath):
        raise SystemExit(f"[erro] Não encontrei {mpath}")
    M = pd.read_csv(mpath, index_col=0)
    pops = list(M.index)

    # --- carregar coordenadas ---
    if not os.path.exists(args.coords):
        raise SystemExit("[erro] Forneça coords em dados/geo_coords.csv (Geo_loc_name,lat,lon).")
    C = pd.read_csv(args.coords)
    C["Geo_loc_name"] = C["Geo_loc_name"].astype(str)

    # alinhar ordem das coordenadas com a matriz Fst
    cc = C.set_index("Geo_loc_name").reindex(pops)
    if cc.isna().any().any():
        faltam = cc[cc.isna().any(axis=1)].index.tolist()
        raise SystemExit(f"[erro] Faltam coordenadas para: {faltam}")

    # quantas populações?
    n_pops = len(pops)
    if n_pops < 2:
        raise SystemExit("[aviso] Menos de 2 populações — análise não aplicável.")

    # --- matriz geográfica ---
    G = np.zeros((n_pops, n_pops), dtype=float)
    for i,p1 in enumerate(pops):
        for j,p2 in enumerate(pops):
            if i<j:
                G[i,j]=G[j,i]=haversine_km(cc.loc[p1,"lat"], cc.loc[p1,"lon"],
                                           cc.loc[p2,"lat"], cc.loc[p2,"lon"])

    # --- preparar matriz F_ST (substituir NaN por média para gráficos; Mantel usa apenas quando >=3) ---
    A = M.values.astype(float)
    np.fill_diagonal(A, 0.0)
    finite = A[np.isfinite(A) & (~np.eye(n_pops, dtype=bool))]
    fill = np.nanmean(finite) if finite.size>0 else 0.5
    A[~np.isfinite(A)] = fill

    # --- montar pontos e rótulos para scatter (apenas entradas definidas em M) ---
    xs, ys, labels = [], [], []
    for i in range(n_pops):
        for j in range(i+1, n_pops):
            if np.isfinite(M.values[i,j]):
                xs.append(G[i,j]); ys.append(M.values[i,j])
                p1, p2 = pops[i], pops[j]
                if args.abbr:
                    p1, p2 = abbr3(p1), abbr3(p2)
                labels.append(f"{p1}–{p2}")
    xs, ys = np.array(xs), np.array(ys)

    # --- caso com 2 populações: sem Mantel, apenas relatório do par ---
    if n_pops == 2:
        out_png = os.path.join(args.species_dir, "mantel_scatter.png")
        plt.figure(figsize=(6.5,5))
        if xs.size == 1 and ys.size == 1:
            plt.scatter(xs, ys, s=60)
            if args.label_points and labels:
                # deslocamento pequeno para evitar sobreposição com o marcador
                plt.text(xs[0]*1.001, ys[0]*1.001, f" {labels[0]}", fontsize=9, va='bottom', ha='left')
            plt.xlabel("Distância geográfica (km)")
            plt.ylabel("F_ST (Hudson)")
            ttl = f"2 populações — distância={xs[0]:.1f} km, F_ST={ys[0]:.3f}"
            plt.title(ttl)
        else:
            plt.text(0.5,0.5,"Sem par F_ST válido", ha='center', va='center')
        plt.tight_layout(); plt.savefig(out_png, dpi=300); plt.close()

        # resumo em txt
        with open(os.path.join(args.species_dir,"mantel_result.txt"), "w", encoding="utf-8") as f:
            f.write("Mantel test: não aplicável (são necessárias >= 3 populações).\n")
            if xs.size == 1 and ys.size == 1:
                f.write(f"Par único: distância = {xs[0]:.3f} km, F_ST = {ys[0]:.6f}\n")
                f.write("Sugestão: incluir >=1 população adicional para testar IBD.\n")
        print("[ok] 2 populações: Mantel não aplicável. Scatter e resumo salvos.")
        return

    # --- caso com >=3 populações: rodar Mantel ---
    try:
        from skbio.stats.distance import mantel
    except Exception as e:
        raise SystemExit("[erro] scikit-bio não encontrado. Instale com: pip install scikit-bio")

    r, p, _ = mantel(A, G, method='pearson', permutations=999, alternative='two-sided')

    # scatter com reta de tendência + rótulos
    out_png = os.path.join(args.species_dir, "mantel_scatter.png")
    plt.figure(figsize=(6.5,5))
    if xs.size >= 2:
        plt.scatter(xs, ys, s=40)
        # regressão linear simples
        b1 = np.cov(xs, ys, ddof=1)[0,1] / np.var(xs, ddof=1)
        b0 = np.mean(ys) - b1*np.mean(xs)
        xline = np.linspace(xs.min(), xs.max(), 100)
        yline = b0 + b1*xline
        plt.plot(xline, yline)
        # rótulos (se pedido)
        if args.label_points and labels:
            for x, y, lab in zip(xs, ys, labels):
                plt.text(x*1.001, y*1.001, f" {lab}", fontsize=8, va='bottom', ha='left')
    plt.xlabel("Distância geográfica (km)")
    plt.ylabel("F_ST (Hudson)")
    plt.title(f"Mantel: r={r:.3f}, p={p:.3f}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300); plt.close()

    with open(os.path.join(args.species_dir,"mantel_result.txt"), "w", encoding="utf-8") as f:
        f.write("Mantel test (Pearson, 999 perms)\n")
        f.write(f"r = {r:.6f}\n")
        f.write(f"p = {p:.6f}\n")

    print(f"[ok] Mantel r={r:.3f}, p={p:.3f}")
    print("[ok] Arquivos salvos em:", args.species_dir)

if __name__ == "__main__":
    main()
