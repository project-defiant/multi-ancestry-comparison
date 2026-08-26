#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
gwas_path <- args[1]
ld_path <- args[2]
sample_size <- as.numeric(args[3])
ancestry <- args[4]
output_path <- args[5]

gwas <- read.delim(gwas_path, sep = "\t", stringsAsFactors = FALSE)
ld <- read.delim(ld_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(all(gwas$snp == colnames(ld)))

R <- as.matrix(ld)
rownames(R) <- colnames(R)
z <- gwas$z

fit <- susieR::susie_rss(z = z, R = R, n = sample_size, L = 1)
cs <- susieR::susie_get_cs(fit, coverage = 0.95)
pip <- susieR::susie_get_pip(fit)

cs_id <- rep(0L, length(z))
if (length(cs$cs) > 0) {
  cs_id[cs$cs[[1]]] <- 1L
}

out <- data.frame(
  snp = gwas$snp,
  pos = gwas$pos,
  z = z,
  pip = pip,
  cs_id = cs_id
)

write.table(out, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

cs_size <- sum(cs_id == 1)
cat(sprintf("[%s] susieR credible set size: %d\n", ancestry, cs_size))
if (cs_size > 0) {
  cat(sprintf("[%s] CS SNPs: %s\n", ancestry, paste(out$snp[cs_id == 1], collapse = ", ")))
}
