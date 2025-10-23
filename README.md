
# Projeto: Filogeografia de *Nasua* (CytB) — pipeline completo

Este README documenta **todo o fluxo** que construímos — do preparo dos dados até os resultados de estrutura, diversidade e geografia — usando **Python**, **R** e algumas checagens no **MEGA‑X**. Os exemplos de comando estão prontos para **Windows (CMD)**.

---

## 0) Estrutura de pastas do projeto

```
nasua_filogeo/
├─ dados/
│  ├─ Nasua_F_ST_aligned.fasta
│  ├─ Nasua_CYTB_metadata_com_localizacao.csv
│  └─ geo_coords.csv                 # (opcional; ver exemplo abaixo)
├─ resultados/
│  ├─ Nasua_nasua/
│  └─ Nasua_narica/
├─ scripts/
│  ├─ contar_populacoes.py           # (A) contagens por população
│  ├─ fst_nasua.py                   # F_ST Hudson por espécie
│  ├─ plot_fst_panels.py             # painéis gráficos além do heatmap
│  ├─ haplo_network_and_map.py       # (B) rede de haplótipos + mapa
│  ├─ amova_phi_st.R                 # (C) AMOVA + Φ_ST
│  ├─ mantel_fst_geo.py              # (C) Mantel + labels de pares
│  └─ diversity_stats.py             # (D) π, Hd, S, Tajima D, Dxy
└─ .venv/ (opcional, Python virtualenv)
```

---

## 1) Dados de entrada

- **FASTA alinhado** (CytB) com cabeçalhos tipo `>Accession_Especie`  
  Ex.: `>MK144305.1_Nasua_nasua` ou `>MK135608.1_Nasua_narica`
- **CSV de metadados** com o cabeçalho:
  ```csv
  Accession,Organism,Country,Geo_loc_name,Collector,Collection_date,Locality,Sequence_length,Definition
  ```
- (Opcional) **Coordenadas** de cada população (centroide por país/localidade) em `dados/geo_coords.csv`:
  ```csv
  Geo_loc_name,lat,lon
  Brazil,-14.2350,-51.9253
  Argentina,-38.4161,-63.6167
  Mexico,23.6345,-102.5528
  Guatemala,15.7835,-90.2308
  Belize,17.1899,-88.4977
  Panama,8.5380,-80.7821
  USA,37.0902,-95.7129
  Costa Rica,9.7489,-83.7534
  ```

> **Importante:** O FASTA **precisa estar alinhado** (mesmo comprimento para todas as sequências).


---

## 2) Ambiente (Python e R)

### Python (recomendado usar `.venv`)
```cmd
python -m venv .venv
.\.venv\Scripts\activate
pip install biopython pandas numpy scikit-allel matplotlib scipy scikit-learn networkx scikit-bio
```

### R com `renv` (reprodutibilidade)
No R, dentro da raiz do projeto:
```r
renv::activate()
install.packages(c("ape","pegas","readr","dplyr","stringr"))
renv::snapshot()
```
Para rodar scripts R:
```cmd
Rscript scripts\amova_phi_st.R --species "Nasua narica"
Rscript scripts\amova_phi_st.R --species "Nasua nasua"
```

---

## 3) (A) Contar indivíduos por população (a partir do FASTA presente)

**Script:** `scripts/contar_populacoes.py`  
**Uso:**
```cmd
python scripts\contar_populacoes.py ^
  --fasta dados\Nasua_F_ST_aligned.fasta ^
  --metadata dados\Nasua_CYTB_metadata_com_localizacao.csv ^
  --pop-field Geo_loc_name ^
  --min-n 3 ^
  --outdir resultados
```

**Saídas:**
- `resultados\counts_by_geoloc.csv` — contagem por população
- `resultados\counts_by_species_geoloc.csv` — contagem por espécie × população

> Exemplo de distribuição (CytB) utilizada:
> - *Nasua narica*: Belize (2), Costa Rica (1), Guatemala (20), Mexico (33), Panama (12), USA (15)  
> - *Nasua nasua*: Argentina (7), Brazil (3)

---

## 4) F_ST (Hudson) por espécie

**Script:** `scripts/fst_nasua.py`  
Calcula **Hudson’s F_ST** par‑a‑par para *N. nasua* e *N. narica* separadamente, definindo populações via `--pop-field` (ex.: `Geo_loc_name`).

**Uso:**
```cmd
python scripts\fst_nasua.py ^
  --fasta dados\Nasua_F_ST_aligned.fasta ^
  --metadata dados\Nasua_CYTB_metadata_com_localizacao.csv ^
  --pop-field Geo_loc_name ^
  --min-n 3 ^
  --resultados resultados
```

**Saídas (por espécie em `resultados/<espécie>/`):**
- `pop_sizes.csv` — tamanho amostral por população
- `fst_pairwise_hudson.csv` — matriz F_ST
- `fst_pair_details.csv` — somatórios (num, den) e nº de sítios
- `fst_pairwise_hudson_heatmap.png` — heatmap

---

## 5) Painéis gráficos além do heatmap (MDS, UPGMA, MST)

**Script:** `scripts/plot_fst_panels.py`  
Gera:
- **MDS/PCoA** (mapa 2D das populações)  
- **Dendrograma UPGMA** (clustering sobre F_ST)  
- **MST** (árvore mínima com pesos = F_ST)  
- (opcional) **heatmap** reprocesado

**Uso:**
```cmd
python scripts\plot_fst_panels.py --species-dir resultados\Nasua_nasua --make-heatmap
python scripts\plot_fst_panels.py --species-dir resultados\Nasua_narica --make-heatmap
```

**Saídas (por espécie):**
- `fst_mds.png`, `fst_dendrogram_upgma.png`, `fst_mst.png`, `fst_heatmap_redo.png`

---

## 6) (B) Rede de haplótipos (MST) + mapa

**Script:** `scripts/haplo_network_and_map.py`  
Agrupa sequências idênticas (haplótipos), constrói **MST** entre haplótipos (diferenças de Hamming) e plota **rede** (fatias por composição de população). Se houver `geo_coords.csv`, plota também **mapa** com tortas por haplótipo.

**Uso:**
```cmd
python scripts\haplo_network_and_map.py --species "Nasua narica" --coords dados\geo_coords.csv
python scripts\haplo_network_and_map.py --species "Nasua nasua"  --coords dados\geo_coords.csv
```

**Saídas (por espécie):**
- `haplotypes.csv` — sequência do haplótipo e tamanho
- `haplotype_assignment.csv` — indivíduo → haplótipo
- `haplo_network.png` — rede de haplótipos
- `haplo_map.png` — mapa por população (se coords disponíveis)

---

## 7) (C) AMOVA e Φ_ST (R, pacote **pegas**)

**Script:** `scripts/amova_phi_st.R` (versão robusta)  
- AMOVA global (`amova_summary.txt`) e **Φ_ST** global  
- **Φ_ST par‑a‑par** (pula pares com **n<2** por população ou sem variação)

**Uso:**
```cmd
Rscript scripts\amova_phi_st.R --species "Nasua narica"
Rscript scripts\amova_phi_st.R --species "Nasua nasua"
```

**Saídas (por espécie):**
- `amova_summary.txt`  
- `phi_st_pairwise.csv`  
- `pop_sizes.csv`

> Dica: Para pops muito pequenas (p.ex., n=1), Φ_ST par‑a‑par sai **NA** (corretamente).

---

## 8) (C) Mantel (F_ST × distância geográfica) — com rótulos

**Script:** `scripts/mantel_fst_geo.py` (com flags de rótulo)  
- Para **≥3 populações**: executa **Mantel** (Pearson, 999 permutações) e plota regressão.  
- Para **2 populações**: reporta o par único (distância e F_ST) e gera o gráfico mesmo assim.  
- Opcional: rótulos nos pontos e abreviações de país.

**Uso:**
```cmd
python scripts\mantel_fst_geo.py --species-dir resultados\Nasua_nasua --label-points
python scripts\mantel_fst_geo.py --species-dir resultados\Nasua_narica --label-points --abbr
```

**Saídas (por espécie):**
- `mantel_scatter.png`  
- `mantel_result.txt`

---

## 9) (D) Métricas de diversidade e divergência (π, Hd, S, Tajima’s D, Dxy)

**Script:** `scripts/diversity_stats.py` (com Tajima’s D protegido)  
- **Por população**: n, nº de haplótipos, **Hd**, **π**, **S**, **Tajima’s D**  
- **Entre populações**: **Dxy** e K médio  
- Gráficos: barras de **π** e **Hd**

**Uso:**
```cmd
python scripts\diversity_stats.py --species both
```

**Saídas (por espécie):**
- `diversity_by_population.csv`  
- `dxy_pairwise.csv`  
- `pi_by_population.png`, `Hd_by_population.png`

> Em populações com **S=0** ou **n muito baixo**, Tajima’s D é **NA** (indefinido), como esperado.

---

## 10) Checagens no MEGA‑X (qualidade/curadoria)

- **Reverse Complement** de sequências invertidas: Alignment Explorer → `Sequence → Reverse Complement` (DNA).  
- **Traduzir** (código genético **mitocondrial de vertebrados**) para detectar **stops** e excluir **numts/contaminação**.  
- **BLAST** direto do MEGA para confirmar espécie/gene.  
- **Aparar** alinhamento para a região **comum** a todas as amostras quando necessário.

---

## 11) Boas práticas e interpretação

- Rodar **análises por marcador** e **por espécie** (aqui: CytB).  
- Evitar comparar **valores absolutos** de F_ST entre marcadores diferentes.  
- **n pequenos** (ex.: *N. nasua* Brasil = 3) → resultados **exploratórios**; tente ampliar amostragem.  
- Reportar **π** e **Hd** junto com F_ST/Φ_ST e redes de haplótipos melhora muito a leitura filogeográfica.

---

## 12) Problemas comuns & soluções rápidas

- **Mantel exige ≥3 populações** → com 2 pops, o script já faz o fallback (sem Mantel, só gráfico e resumo).  
- **Φ_ST par‑a‑par falha** com n<2 ou sem variação → nossa versão retorna **NA**.  
- **Tajima’s D warning** (sqrt inválido) → corrigido no script (retorna **NA**).  
- **FASTA não alinhado** → alinhar (MAFFT) e manter apenas região sobreposta.

---

## 13) Referências (insights)

- Artigos/Tese anexados ao projeto (consultados para orientar escolhas analíticas):  
  - `nasua_article_1.pdf`  
  - `nasua_article_2.pdf`  
  - `tese_nasua.pdf`

---

## 14) Resumo dos comandos (Windows)

```cmd
:: (A) contagens por população
python scripts\contar_populacoes.py --fasta dados\Nasua_F_ST_aligned.fasta --metadata dados\Nasua_CYTB_metadata_com_localizacao.csv --pop-field Geo_loc_name --min-n 3 --outdir resultados

:: F_ST por espécie
python scripts\fst_nasua.py --fasta dados\Nasua_F_ST_aligned.fasta --metadata dados\Nasua_CYTB_metadata_com_localizacao.csv --pop-field Geo_loc_name --min-n 3 --resultados resultados

:: painéis gráficos além do heatmap
python scripts\plot_fst_panels.py --species-dir resultados\Nasua_nasua --make-heatmap
python scripts\plot_fst_panels.py --species-dir resultados\Nasua_narica --make-heatmap

:: (B) rede de haplótipos + mapa
python scripts\haplo_network_and_map.py --species "Nasua nasua"  --coords dados\geo_coords.csv
python scripts\haplo_network_and_map.py --species "Nasua narica" --coords dados\geo_coords.csv

:: (C) AMOVA / Φ_ST (R + renv)
Rscript scripts\amova_phi_st.R --species "Nasua nasua"
Rscript scripts\amova_phi_st.R --species "Nasua narica"

:: (C) Mantel (com rótulos)
python scripts\mantel_fst_geo.py --species-dir resultados\Nasua_nasua --label-points
python scripts\mantel_fst_geo.py --species-dir resultados\Nasua_narica --label-points --abbr

:: (D) diversidade e divergência
python scripts\diversity_stats.py --species both
```

---

**Contato/Notas:**  
- Scripts preparados para *Nasua nasua* e *Nasua narica* com **CytB**; é fácil extender para outros marcadores mantendo o paradigma.
- Para publicação, podemos gerar um **relatório PDF/Quarto** integrando tabelas + figuras automaticamente.
