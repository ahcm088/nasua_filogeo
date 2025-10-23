# amova_phi_st.R  —  AMOVA e Φ_ST (pegas) por espécie a partir de FASTA alinhado
suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(readr)
  library(dplyr)
  library(stringr)
})

# -------------------- parsing de argumentos --------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}

fasta_path   <- get_arg("--fasta",   "dados/Nasua_F_ST_aligned.fasta")
metadata_csv <- get_arg("--metadata","dados/Nasua_CYTB_metadata_com_localizacao.csv")
species      <- get_arg("--species", "Nasua narica")     # troque para "Nasua nasua" quando rodar a outra espécie
pop_field    <- get_arg("--pop-field", "Geo_loc_name")
out_dir      <- file.path("resultados", gsub(" ", "_", species))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("[info] FASTA:   ", fasta_path)
message("[info] Meta:    ", metadata_csv)
message("[info] Espécie: ", species)
message("[info] PopField:", pop_field)
message("[info] OutDir:  ", out_dir)

# -------------------- leitura e preparo --------------------
# FASTA alinhado (lista DNAbin)
aln <- read.FASTA(fasta_path)
if (length(aln) == 0) stop("FASTA vazio.")

# Metadados
md  <- read_csv(metadata_csv, show_col_types = FALSE)

# Cabeçalhos -> Accession (ex.: "MK144305.1_Nasua_nasua" -> "MK144305.1")
ids <- names(aln)
acc <- str_replace(ids, "_.*$", "")
org_from_header <- str_replace(ids, "^[^_]+_", "") %>% str_replace_all("_"," ")

df <- tibble(header = ids, Accession = acc, Org_hdr = org_from_header) %>%
      left_join(md, by = "Accession")

# Filtra espécie
df_sp <- df %>% filter(Organism == species)
if (nrow(df_sp) == 0) stop(paste("Sem sequências após merge para:", species))

# Subconjunto do alinhamento na ordem de df_sp
aln_sp <- aln[df_sp$header]     # já é DNAbin
dna    <- aln_sp

# Populações (fator)
if (!pop_field %in% names(df_sp)) {
  stop(sprintf("Campo '%s' não encontrado no metadata.", pop_field))
}
pops <- df_sp[[pop_field]] |> as.character() |> str_trim()
pops[pops %in% c("", "NA")] <- NA
pops <- as.factor(pops)

# Salvar tamanhos amostrais por população (ignora NA)
tab_p <- table(pops, useNA = "no")
pop_sizes <- data.frame(
  populacao = names(tab_p),
  n = as.integer(tab_p),
  row.names = NULL
)
write.csv(pop_sizes, file.path(out_dir, "pop_sizes.csv"), row.names = FALSE)

# -------------------- distâncias e AMOVA global --------------------
# Distância nucleotídica simples (pairwise deletion)
D <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = FALSE)

# AMOVA global: dist ~ pops
amv <- amova(D ~ pops)
phi_st_global <- amv$phi[1, "Phi"]

sink(file.path(out_dir, "amova_summary.txt"))
cat("=== AMOVA (pegas) ===\n"); print(amv)
cat("\nPhiST global:", phi_st_global, "\n")
sink()

message(sprintf("[ok] AMOVA concluída. PhiST global = %.4f", phi_st_global))

# -------------------- Φ_ST par-a-par --------------------
lev <- levels(pops)
# manter apenas níveis realmente presentes (descarta nível NA)
lev <- lev[lev %in% pop_sizes$populacao]

pairs <- combn(lev, 2, simplify = FALSE)

res <- lapply(pairs, function(pp){
  idx <- which(pops %in% pp)
  if (length(idx) < 2) return(NULL)

  sub_pops <- droplevels(pops[idx])
  n_tab <- table(sub_pops)

  # cada população do par precisa ter pelo menos 2 amostras para estimativa estável
  if (any(n_tab < 2)) {
    return(data.frame(pop1 = pp[1], pop2 = pp[2], PhiST = NA_real_,
                      n1 = unname(n_tab[1]),
                      n2 = ifelse(length(n_tab) > 1, unname(n_tab[2]), NA_integer_)))
  }

  # submatriz de distâncias para os índices selecionados
  Dm <- as.matrix(D)
  Dp <- as.dist(Dm[idx, idx])

  # se não há variação (todas as distâncias 0), Φ_ST é indefinido
  if (all(Dp == 0)) {
    return(data.frame(pop1 = pp[1], pop2 = pp[2], PhiST = NA_real_,
                      n1 = unname(n_tab[1]), n2 = unname(n_tab[2])))
  }

  # AMOVA para o par com proteção a casos degenerados
  amv2 <- tryCatch(pegas::amova(Dp ~ sub_pops), error = function(e) NULL)

  if (is.null(amv2) || is.null(amv2$phi) || nrow(amv2$phi) < 1 || is.na(amv2$phi[1, "Phi"])) {
    return(data.frame(pop1 = pp[1], pop2 = pp[2], PhiST = NA_real_,
                      n1 = unname(n_tab[1]), n2 = unname(n_tab[2])))
  }

  data.frame(pop1 = pp[1], pop2 = pp[2],
             PhiST = amv2$phi[1, "Phi"],
             n1 = unname(n_tab[1]), n2 = unname(n_tab[2]))
})

res <- do.call(rbind, res)
if (!is.null(res) && nrow(res) > 0) {
  write.csv(res, file.path(out_dir, "phi_st_pairwise.csv"), row.names = FALSE)
  message("[ok] PhiST par-a-par salvo em ", file.path(out_dir, "phi_st_pairwise.csv"))
} else {
  message("[aviso] Não foi possível calcular pares informativos (amostra insuficiente).")
}

message("[feito] Saídas em: ", out_dir)
