#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS

def load_fst_matrix(species_dir):
    mpath = os.path.join(species_dir, "fst_pairwise_hudson.csv")
    if not os.path.exists(mpath):
        raise SystemExit(f"[erro] arquivo não encontrado: {mpath}")
    M = pd.read_csv(mpath, index_col=0)
    # Garantir simetria e diagonal 0
    pops = list(M.index)
    M = M.loc[pops, pops]  # reordenar igual ao índice
    np.fill_diagonal(M.values, 0.0)
    # Simetrizar
    M = (M + M.T) / 2.0
    return M

def load_pop_sizes(species_dir):
    ppath = os.path.join(species_dir, "pop_sizes.csv")
    if os.path.exists(ppath):
        df = pd.read_csv(ppath)
        sizes = dict(zip(df["populacao"].astype(str), df["n"].astype(int)))
        return sizes
    return {}

def impute_nan_distances(M):
    """Imputa NaN off-diagonal com uma penalidade (maior Fst observado + 0.05)."""
    A = M.values.astype(float).copy()
    # diagonal zero
    np.fill_diagonal(A, 0.0)
    off_diag = A[~np.eye(A.shape[0], dtype=bool)]
    finite_vals = off_diag[np.isfinite(off_diag)]
    if finite_vals.size == 0:
        # tudo NaN (caso extremo)
        A[~np.eye(A.shape[0], dtype=bool)] = 1.0
        return pd.DataFrame(A, index=M.index, columns=M.columns)
    penalty = np.nanmax(finite_vals) + 0.05
    # substitui NaN fora da diagonal pela penalidade
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i != j and not np.isfinite(A[i, j]):
                A[i, j] = penalty
    # força simetria
    A = (A + A.T) / 2.0
    return pd.DataFrame(A, index=M.index, columns=M.columns)

def plot_heatmap(M, out_png, title="Hudson F_ST (pairwise)"):
    plt.figure(figsize=(7, 6))
    im = plt.imshow(M.values.astype(float), interpolation='nearest')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    pops = list(M.index)
    plt.xticks(range(len(pops)), pops, rotation=45, ha='right')
    plt.yticks(range(len(pops)), pops)
    plt.title(title)
    # Anotar
    for i in range(len(pops)):
        for j in range(len(pops)):
            if i == j:
                continue
            val = M.values[i, j]
            if np.isfinite(val):
                plt.text(j, i, f"{val:.3f}", ha='center', va='center', fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_mds(M, pop_sizes, out_png, title="MDS (a partir de F_ST)"):
    # MDS métrico em dissimilaridades
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=0, n_init=8, max_iter=1000)
    X = mds.fit_transform(M.values)
    pops = list(M.index)
    sizes = np.array([pop_sizes.get(p, 3) for p in pops])
    # escala tamanho (opcional)
    s = 80 + 20 * (sizes - sizes.min()) if sizes.size and sizes.max() > sizes.min() else 120

    plt.figure(figsize=(7, 6))
    plt.scatter(X[:, 0], X[:, 1], s=s)
    for i, p in enumerate(pops):
        plt.text(X[i, 0], X[i, 1], f" {p}", va='center', ha='left', fontsize=10)
    plt.title(title)
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_dendrogram_upgma(M, out_png, title="Dendrograma (UPGMA sobre F_ST)"):
    # Condensed vector para clustering
    D = squareform(M.values, checks=False)  # usa somente triângulo superior
    Z = hierarchy.linkage(D, method='average')  # UPGMA ~ average linkage
    plt.figure(figsize=(8, 5))
    hierarchy.dendrogram(Z, labels=list(M.index), orientation='right', leaf_font_size=10)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_mst(M, pop_sizes, out_png, title="Árvore de Mínima Abrangência (MST) por F_ST"):
    pops = list(M.index)
    G = nx.Graph()
    for p in pops:
        G.add_node(p, size=pop_sizes.get(p, 3))

    # adicionar arestas para pares com valor finito
    for i in range(len(pops)):
        for j in range(i+1, len(pops)):
            w = M.values[i, j]
            if np.isfinite(w):
                G.add_edge(pops[i], pops[j], weight=float(w))

    if G.number_of_edges() == 0:
        print("[aviso] Sem arestas finitas para MST.")
        return

    # MST
    T = nx.minimum_spanning_tree(G, weight='weight')
    # layout usando pesos (distâncias) -> converte pra proximidade para spring_layout
    # quanto menor Fst, mais "forte" a mola
    inv_w = {e: 1.0 / (1.0 + d['weight']) for e, d in T.edges.items()}
    pos = nx.spring_layout(T, weight=None, k=None, seed=42)  # layout básico

    plt.figure(figsize=(8, 6))
    # nós
    node_sizes = [100 + 30 * (G.nodes[p]['size']) for p in T.nodes]
    nx.draw_networkx_nodes(T, pos, node_size=node_sizes)
    # arestas (espessura ~ proximidade)
    widths = [3.0 * inv_w[e] for e in T.edges]
    nx.draw_networkx_edges(T, pos, width=widths)
    # rótulos
    nx.draw_networkx_labels(T, pos, font_size=10)
    # rótulos nas arestas com F_ST
    edge_labels = {(u, v): f"{T[u][v]['weight']:.3f}" for u, v in T.edges}
    nx.draw_networkx_edge_labels(T, pos, edge_labels=edge_labels, font_size=8)

    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def main():
    ap = argparse.ArgumentParser(description="Gera visualizações além do heatmap para matrizes de F_ST (Hudson).")
    ap.add_argument("--species-dir", required=True, help="Pasta da espécie (ex.: resultados/Nasua_nasua)")
    ap.add_argument("--make-heatmap", action="store_true", help="Se passado, (re)gera um heatmap.")
    args = ap.parse_args()

    os.makedirs(args.species_dir, exist_ok=True)

    # Carregar dados
    M_raw = load_fst_matrix(args.species_dir)
    pop_sizes = load_pop_sizes(args.species_dir)

    # Imputar NaNs para métodos que não aceitam NaN (MDS/dendrograma/MST)
    M = impute_nan_distances(M_raw)

    # 1) (opcional) Heatmap
    if args.make_heatmap:
        plot_heatmap(M_raw, os.path.join(args.species_dir, "fst_heatmap_redo.png"),
                    title="Hudson F_ST (pairwise)")

    # 2) MDS
    plot_mds(M, pop_sizes, os.path.join(args.species_dir, "fst_mds.png"),
             title="MDS (a partir de F_ST)")

    # 3) Dendrograma (UPGMA)
    if len(M.index) >= 2:
        plot_dendrogram_upgma(M, os.path.join(args.species_dir, "fst_dendrogram_upgma.png"),
                              title="Dendrograma (UPGMA sobre F_ST)")

    # 4) MST
    if len(M.index) >= 2:
        plot_mst(M, pop_sizes, os.path.join(args.species_dir, "fst_mst.png"),
                 title="Árvore de Mínima Abrangência (MST) por F_ST")

    print("[ok] Gráficos salvos em:", args.species_dir)

if __name__ == "__main__":
    main()
