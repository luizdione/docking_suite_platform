#!/usr/bin/env Rscript
# =============================================================================
# run_dep_limma.R - enriquecimento diferencial AP-MS em nivel de proteina.
#
# Usa DEP (Zhang et al., Nat Protoc 2018), que envolve limma com filtragem,
# normalizacao vsn e imputacao MNAR (MinProb) - o fluxo de referencia para
# matrizes de proteinGroups/pg_matrix/combined_protein. NAO reimplementa nada:
# apenas orquestra DEP/limma.
#
# Uso:
#   Rscript run_dep_limma.R data.tsv design.tsv out.tsv <test> <control> [alpha] [lfc]
#
# data.tsv   : colunas ID, name e uma coluna de intensidade por amostra (label)
# design.tsv : colunas label, condition, replicate
# out.tsv    : saida (get_results do DEP)
# =============================================================================
suppressMessages({
  library(DEP)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Uso: run_dep_limma.R data.tsv design.tsv out.tsv <test> <control> [alpha] [lfc]")
}
data_path   <- args[1]
design_path <- args[2]
out_path    <- args[3]
test_cond   <- args[4]
ctrl_cond   <- args[5]
alpha       <- ifelse(length(args) >= 6, as.numeric(args[6]), 0.05)
lfc         <- ifelse(length(args) >= 7, as.numeric(args[7]), 1.0)

data <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
design <- read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)

# nomes unicos exigidos pelo DEP
data_unique <- make_unique(data, "name", "ID", delim = ";")

# colunas de intensidade = labels do design
lfq_columns <- match(design$label, colnames(data_unique))
if (any(is.na(lfq_columns))) {
  stop("Labels do design ausentes em data.tsv: ",
       paste(design$label[is.na(lfq_columns)], collapse = ", "))
}

se <- make_se(data_unique, lfq_columns, design)

# filtragem: proteina presente (sem faltantes) em ao menos uma condicao
data_filt <- filter_missval(se, thr = 0)
# normalizacao vsn (estabiliza a variancia)
data_norm <- normalize_vsn(data_filt)
# imputacao MNAR (left-censored) apropriada para AP-MS
data_imp <- tryCatch(
  impute(data_norm, fun = "MinProb", q = 0.01),
  error = function(e) { message("impute MinProb falhou, usando 'man': ", conditionMessage(e));
                        impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) }
)

# contraste manual test_vs_control (nomes = condition do design)
contrast <- paste0(test_cond, "_vs_", ctrl_cond)
data_diff <- test_diff(data_imp, type = "manual", test = contrast)
dep <- add_rejections(data_diff, alpha = alpha, lfc = lfc)

res <- get_results(dep)
write.table(res, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("run_dep_limma: %d proteinas, contraste %s, alpha=%.3g lfc=%.3g -> %s\n",
            nrow(res), contrast, alpha, lfc, out_path))
