
# INSTRUÇÕES — Interpretação das saídas do pipeline de filogeografia (*Nasua*, CytB)

Este documento explica **a técnica por trás** de cada etapa e dá **instruções de interpretação** para **todas as saídas** listadas no `README.md` — do censo por população às análises de estrutura (F_ST/Φ_ST/Mantel), redes de haplótipos, e métricas de diversidade (π, Hd, S, Tajima’s D, Dxy).

> Use este guia lado a lado com o `README.md` para saber **o que cada arquivo/figura significa** e **como reportar** os resultados no texto do estudo.


---

## 1) Dados e pré‑processamento

### 1.1 FASTA alinhado (`dados/Nasua_F_ST_aligned.fasta`)
- **O que é**: sequências CytB **já alinhadas** e de **mesmo comprimento**.
- **Por que importa**: quase todas as estatísticas (F_ST, π, Dxy, redes) assumem **comparação posição‑a‑posição**. Misalignments introduzem diferenças falsas.
- **Checklist de qualidade**:
  - Verifique no MEGA‑X: sem *stop codons* (tradução com **código mitocondrial de vertebrados**), sem regiões inversas (se houver, aplique *Reverse Complement*), sem excesso de `N` ou *gaps*.
  - **Apare** regiões não sobrepostas: mantenha o *core* comum a todas as amostras.

### 1.2 Metadados (`dados/Nasua_CYTB_metadata_com_localizacao.csv`)
- **Chaves**: `Accession` (casa com o FASTA), `Organism` (espécie), `Geo_loc_name` (população/país).
- **Boas práticas**: padronize `Geo_loc_name` (ex.: “Mexico” ≠ “México”), mantenha uma linha por indivíduo.


---

## 2) (A) Contagem por população

### `counts_by_geoloc.csv` e `counts_by_species_geoloc.csv`
- **O que é**: censo de amostras por **população** e por **espécie × população**.
- **Como usar**:
  - Avalie **poder amostral**. Estimativas com **n < 5** por população são ruidosas.
  - Decida **quais comparações são válidas** (ex.: pares com n≥2 em cada lado para Φ_ST).


---

## 3) F_ST de Hudson (estrutura genética)

### 3.1 `fst_pairwise_hudson.csv` (matriz) e `fst_pair_details.csv` (componentes)
- **Técnica**: F_ST de **Hudson** baseado em **π_w** (dentro) e **π_b** (entre). É robusto para sequências (mtDNA) e tamanhos desiguais.
- **Interpretação**:
  - **Valores**: 0 ⇒ populações indistintas; >0 ⇒ estrutura; aproximar de 1 ⇒ altamente distintas.
  - **Atenção**: F_ST depende da **diversidade interna**. Populações com π muito baixo tendem a inflar F_ST.
  - **Use junto com**: **Φ_ST** (AMOVA), **Dxy** e **rede de haplótipos**.

### 3.2 `fst_pairwise_hudson_heatmap.png`
- **O que olhar**: blocos/quadrantes mais claros (F_ST alto) indicam **barreiras**; gradientes sugerem **IBD** (*isolation by distance*).
- **Boas práticas**: cite **n por população** (arquivo `pop_sizes.csv`) ao discutir contrastes fortes/ fracos.


---

## 4) Painéis além do heatmap (estrutura e relações)

### 4.1 `fst_mds.png` (MDS/PCoA a partir de F_ST)
- **Técnica**: projeta a matriz de dissimilaridade (F_ST) em 2D para visualizar **afinidades**.
- **Leitura**: clusters próximos ⇒ populações geneticamente similares; eixos **não** têm significado direto biológico.
- **Cuidado**: MDS requer matriz **completa**; no script imputamos valores faltantes com penalidade — descreva no método.

### 4.2 `fst_dendrogram_upgma.png` (UPGMA sobre F_ST)
- **Técnica**: *average linkage* (UPGMA) na matriz F_ST.
- **Leitura**: ramos curtos ⇒ proximidade; cortes naturais sugerem **grupos**.
- **Limitação**: dendrograma de **distância**, não é árvore filogenética de genes/espécies.

### 4.3 `fst_mst.png` (árvore de mínima abrangência por F_ST)
- **Técnica**: MST com pesos = F_ST (arestas mínimas que conectam todas as pops).
- **Leitura**: caminhos/“pontes” de menor F_ST entre grupos; espessura/rotulagem das arestas auxiliam na leitura de **barreiras**.


---

## 5) (B) Rede de haplótipos e mapa

### 5.1 `haplotypes.csv` e `haplotype_assignment.csv`
- **Técnica**: colapsa sequências **idênticas** em **haplótipos**; contabiliza `Size` e mapeia indivíduo→haplótipo.
- **Uso**: reporte **nº de haplótipos** por população e haplótipos **compartilhados** entre populações.

### 5.2 `haplo_network.png` (MST de haplótipos, fatias por população)
- **Leitura**: nós grandes = haplótipos comuns/centrais; arestas longas = maior nº de mutações; “estrelas” sugerem **expansão recente**.
- **Perguntas-chave**: há **haplótipos exclusivos** por país? Há **linhagens** segregadas?

### 5.3 `haplo_map.png` (tortas por haplótipo no espaço)
- **Leitura**: distribuição espacial de haplótipos; tortas mistas sugerem **fluxo**/contato; tortas “puras” sugerem **isolamento**.
- **Cuidado**: granulação espacial depende do que entrou em `Geo_loc_name` e de `geo_coords.csv` (centroides).


---

## 6) (C) AMOVA e Φ_ST

### 6.1 `amova_summary.txt`
- **Técnica**: AMOVA (*Analysis of Molecular Variance*) com **distância de sequência** (pegas) → particiona variância **entre** e **dentro** de populações.
- **Saídas**: componentes de variância e **Φ_ST global** (análogo a F_ST, mas pondera **distâncias** entre haplótipos).
- **Interpretação**: Φ_ST alto e significativo ⇒ **estrutura** consistente; compare com F_ST.

### 6.2 `phi_st_pairwise.csv`
- **Técnica**: AMOVA 2‑a‑2; ignoramos pares com **n < 2** por população ou **sem variação** (Φ_ST = NA).
- **Uso**: reforça **pares críticos** (barreiras) e serve de “replicação” do padrão visto no F_ST.


---

## 7) (C) Mantel: F_ST × distância geográfica

### 7.1 `mantel_scatter.png` e `mantel_result.txt`
- **Técnica**: correlação entre **distância genética** (F_ST) e **distância geográfica** (km, *haversine*). Permutações (999) para *p‑valor*.
- **Leitura**:
  - **r > 0 e p < 0.05**: sinal de **IBD** (quanto mais longe, mais diferente).
  - **r ≈ 0**: estrutura **não** explicada por distância (talvez barreiras, história).
- **Observações**:
  - Requer **≥ 3 populações**. Com 2 pops, o script faz *fallback* (apenas ponto único + resumo).  
  - Rótulos de pares ajudam a identificar **outliers** (pontos fora da tendência).


---

## 8) (D) Diversidade e divergência

### 8.1 `diversity_by_population.csv` (por população)
- **Métricas**:
  - **n**: tamanho amostral.
  - **n_hap**: número de haplótipos.
  - **Hd** (*haplotype diversity*): probabilidade de dois indivíduos aleatórios terem **haplótipos diferentes**. Próximo de 1 ⇒ alta diversidade de haplótipos.
  - **π** (*nucleotide diversity*): diferenças **por sítio** dentro da população.
  - **S**: nº de **sítios segregantes** (variáveis).
  - **Tajima’s D**: compara π e θ_W; **D < 0** sugere **expansão recente** ou seleção purificadora; **D > 0** sugere **estrangulamento**/estrutura/seleção balanceadora. **NA** quando indefinido (p.ex., S = 0 ou n muito baixo).
- **Interpretação**:
  - Reporte **Hd e π** com **IC/erro** quando possível (pode-se *bootstrap* no futuro).
  - Em populações com **n pequeno** (ex.: 3), trate como **exploratório**.

### 8.2 `dxy_pairwise.csv` (entre populações)
- **Dxy**: diferenças **por sítio** **entre** populações (absoluto, independe da variação interna).
- **Uso**: complementa F_ST; pares com **Dxy alto** mas **F_ST moderado** podem indicar **diversidade interna alta** em ambos os lados.

### 8.3 `pi_by_population.png` e `Hd_by_population.png`
- **Leitura**: barras por população; combine com `pop_sizes.csv` para evitar sobre‑interpretação com **n pequeno**.


---

## 9) Integração e reporte

- **Conte a história** usando **convergência de evidências**:
  1) **Rede de haplótipos + mapa** (linhagens/compartilhamento)  
  2) **F_ST/Φ_ST** (estrutura)  
  3) **Mantel** (IBD vs barreiras)  
  4) **π/Hd/Tajima’s D** (história demográfica)  
  5) **Dxy** (divergência absoluta).
- **Gráficos‑chave para artigo**: rede de haplótipos no mapa; heatmap + MDS; barras de π/Hd; Mantel com rótulos.
- **Limitações**: CytB é **mtDNA** (herança materna); resultados podem não refletir o genoma nuclear. Declarem isso e evitem extrapolações indevidas.


---

## 10) Boas práticas e controles

- **Curadoria**: filtre sequências com muitos `N`/gaps, *stops* ao traduzir CytB, e Numts suspeitos (BLAST/MEGA‑X).
- **Aparar região comum** antes de calcular estatísticas.
- **Balanceamento**: quando comparar países, considere **rarefação** de *N. narica* para n semelhante entre países.
- **Reprodutibilidade**: use `renv` para R e `.venv` para Python; *snapshot/lock* sempre que mudar dependências.


---

## 11) Caminhos para avançar (opcional)

- **SAMOVA** (Arlequin): grupos espaciais que maximizam Φ_ST (K=2..5).
- **Rede TCS** (PopART) além da MST.
- **BPEC** (R): clustering filogeográfico bayesiano (mtDNA + coordenadas).
- **BEAST**: datas de divergência (clock relaxado, HKY+Γ, partição por códons; prior de taxa de CytB para carnívoros).


---

## 12) Anexos úteis

- `pop_sizes.csv` (por espécie): sempre cite ao interpretar F_ST/Φ_ST/Mantel.
- `README.md`: instruções operacionais e comandos do pipeline.
- Artigos/tese anexos: contextualizam expectativas para *Nasua* (linhagens, distribuição).
