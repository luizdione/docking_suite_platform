#!/usr/bin/env Rscript
# =============================================================================
# run_msstats.R - quantificacao diferencial em nivel de FEATURE com MSstats.
#
# MSstats (Choi et al., Bioinformatics 2014) e o padrao quando ha dados em nivel
# de peptideo/feature (evidence.txt do MaxQuant ou report.tsv do DIA-NN). Use os
# conversores oficiais para gerar o input anotado antes deste script:
#   MaxQuant : MSstats::MaxQtoMSstatsFormat(evidence, annotation, proteinGroups)
#   DIA-NN   : MSstatsConvert::DIANNtoMSstatsFormat(report, annotation)
#
# Este script recebe o CSV ja no formato MSstats e roda dataProcess +
# groupComparison para o contraste test-control. NAO reimplementa o modelo.
#
# Uso:
#   Rscript run_msstats.R msstats_input.csv contrasts.tsv out.tsv
#
# contrasts.tsv : colunas 'test' e 'control' (nomes de Condition), uma linha
# =============================================================================
suppressMessages({
  library(MSstats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: run_msstats.R msstats_input.csv contrasts.tsv out.tsv")
}
input_path     <- args[1]
contrasts_path <- args[2]
out_path       <- args[3]

raw <- read.csv(input_path, stringsAsFactors = FALSE)
contrasts <- read.delim(contrasts_path, stringsAsFactors = FALSE)

processed <- dataProcess(
  raw,
  normalization = "equalizeMedians",
  summaryMethod = "TMP",
  censoredInt = "NA",
  MBimpute = TRUE
)

# matriz de contraste test - control sobre os niveis de Condition
conditions <- levels(as.factor(processed$ProteinLevelData$GROUP))
cmat <- matrix(0, nrow = nrow(contrasts), ncol = length(conditions),
               dimnames = list(NULL, conditions))
row_names <- character(nrow(contrasts))
for (i in seq_len(nrow(contrasts))) {
  cmat[i, contrasts$test[i]]    <- 1
  cmat[i, contrasts$control[i]] <- -1
  row_names[i] <- paste0(contrasts$test[i], "_vs_", contrasts$control[i])
}
rownames(cmat) <- row_names

comp <- groupComparison(contrast.matrix = cmat, data = processed)
res <- comp$ComparisonResult
write.table(res, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("run_msstats: %d linhas de comparacao -> %s\n", nrow(res), out_path))
