

library(tidyverse)
library(data.table)
library(glue)
library(cowplot)
library(ggpointdensity)

if (interactive()) {
  tissue  <- "Liver"
  basepath <- "/project/yangili1/cdai/SpliFi"
  out_prefix <- "/project/yangili1/cdai/splice-pub/smk-plots/gtex-sqtl-enrichment"
  FDR <- 0.1
} else {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 4) {
    stop("Usage: Rscript plot-gtex-sqtl-enrichment.R tissue basepath out_prefix FDR")
  }

  tissue <- args[1]
  basepath <- args[2]
  out_prefix <- args[3]
  FDR <- as.numeric(args[4])
}

print(glue("tissue: {tissue}, basepath: {basepath}, out_prefix: {out_prefix}"))





#------ functions   ------

multiqq <- function(pvalues, flipY = FALSE) {
  # flipgY: if TRUE, flip the y-axis
  library(foreach)
  if (is.null(names(pvalues))) {
    names(pvalues) <- seq_along(pvalues)
  }
  punif <- -log10(runif(max(sapply(pvalues, length))))
  df <- do.call(rbind, foreach(i = seq_len(length(pvalues))) %do% {
    df <- as.data.frame(
      qqplot(
        x = punif[1:length(pvalues[[i]])],
        y = -log10(pvalues[[i]]),
        plot.it = FALSE
      )
    )
    df$group <- names(pvalues)[i]
    df
  })

  if (flipY) {
    df$y <- -df$y
  }

  df$group <- factor(df$group, names(pvalues))
  ggplot(df, aes(x, y, col = group)) +
    geom_point() +
    geom_abline(
      intercept = 0,
      slope = 1
    ) +
    xlab("Expected -log10(p)") +
    ylab("Observed -log10(p)")
}

addLabels <- function(dt) {
  pid_split <- str_split(dt$phenotype_id, ":")
  dt$itype <- map_chr(pid_split, ~ .x[5])
  dt$clu <- map_chr(pid_split, ~ .x[4])
  dt$clu <- str_extract(dt$clu, "clu_[0-9]+")

  dt <- dt[, ctype := paste(sort(unique(itype)), collapse = ","), by = clu][]
  return(dt)
}



readGTExSQTL <- function(tissue, basepath) {
  base1.gtex <- glue("{basepath}/code/results/qtl/noisy/GTEx/")
  base2.gtex <- "/cis_100000/perm"
  suffix <- "addQval.txt.gz"
  folder <- glue("{base1.gtex}/{tissue}/separateNoise{base2.gtex}")
  files <- glue("{folder}/chr{1:22}.{suffix}")
  dt <- map_dfr(files, fread)
}

readGTExNOM <- function(tissue, basepath) {
  base1.gtex <- glue("{basepath}/code/results/eqtl/GTEx")
  paths <- glue("{base1.gtex}/{tissue}/nom/chr{1:22}.txt.gz")
  df <- map_dfr(paths, fread)
  names(df) <- c(
    "pid", "pchr", "pstart", "pend", "pstrand", "nVar", "dist",
    "vid", "vchr", "vstart", "vend", "pval", "r2", "slope", "topflag"
  )
  return(df)
}

readGTExNOM_allSnps <- function(tissue, basepath) {
  base1.gtex <- glue("{basepath}/code/results/eqtl/GTEx")
  paths <- glue("{base1.gtex}/{tissue}/nom-all-snps/chr{1:22}.txt.gz")
  df <- map_dfr(paths, fread)
  names(df) <- c(
    "pid", "pchr", "pstart", "pend", "pstrand", "nVar", "dist",
    "vid", "vchr", "vstart", "vend", "pval", "r2", "slope", "topflag"
  )
  # randomly select 2000 snps
  set.seed(123)
  df <- df[sample(nrow(df), 2000), ]
  return(df)
}

# ------ main ------

permDF <- readGTExSQTL(tissue, basepath) %>%
  addLabels
print("dim of permDF:")
print(dim(permDF))

nomDF <- readGTExNOM(tissue, basepath)
print("dim of nomDF")
print(dim(nomDF))

# genomewide snps, randomly select 2000
nomDF_genome <- readGTExNOM_allSnps(tissue, basepath)
print("dim of nomDF_geonomewide")
print(dim(nomDF_genome))


# determine p-sqTLs and u-sQTLs
sqtl <- permDF[q < FDR & itype %in% c("PR", "UP") & ctype %in% c("PR", "PR,UP")]

# select the best intron per cluster based on the largest effect size
# note u-sQTLs are selected based on teh largest unproductive effect size
sqtl <- sqtl[, rk := frank(-abs(slope), ties.method = "first"), by = .(clu, itype)]
sqtl <- rbind(
  sqtl[ctype == "PR"][itype == 'PR' & rk == 1][order(clu)], # best intron for PR cluster
  sqtl[ctype == "PR,UP"][itype == "UP" & rk == 1][order(clu)] # best intron for PR,UP cluster
) %>%
  .[naturalsort::naturalorder(phenotype_id)]

# give key for joining with eqtl nominal pass
sqtl.v <- sqtl[, best_genotype_id] %>% unique()


# first subset nominal data using potentila keys to join
nomDF <- nomDF[vid %in% sqtl.v]

# join tables
mergeDF <- inner_join(sqtl[, .(phenotype_id, phenotype_chr, phenotype_start, phenotype_end, phenotype_strand, best_genotype_id, best_nom_dist, pval_nom, slope, itype, ctype)],
  nomDF[, .(pid, pchr, pstart, pend, pstrand, vid, dist, pval, slope, topflag)],
  by = c("best_genotype_id" = "vid"),
  suffix = c("_sqtl", "_eqtl"),
  relationship = "many-to-many"
)

# a eQTL (gene) must encompass the sQTL (intron)
# thus, the start and end of the eQTL must be within the start and end of the sQTL
# the strand must be the same, and the chromosome must be the same

mergeDF <- mergeDF[phenotype_chr == pchr & phenotype_strand == pstrand]
mergeDF <- mergeDF[phenotype_start >= pstart & phenotype_end <= pend]


#---- plot qqplot ----

plotDF <- mergeDF[, .(
  gid = best_genotype_id,
  intron_id = phenotype_id,
  gene_id = pid,
  pval_sqtl = pval_nom,
  pval_eqtl = pval,
  slope_sqtl,
  slope_eqtl,
  ctype,
  topflag
)]

print("dim of plot data:")
print(dim(plotDF))

# first plot sQTL_slope > 0
qqplot <- list(
  Productive = plotDF[ctype == "PR" & slope_sqtl > 0, pval_eqtl],
  Unproductive = plotDF[ctype == "PR,UP" & slope_sqtl > 0, pval_eqtl]
) %>%
  multiqq()

# then plot sQTL_slope < 0
qqplotNeg <- list(
  Productive = plotDF[ctype == "PR" & slope_sqtl < 0, pval_eqtl],
  Unproductive = plotDF[ctype == "PR,UP" & slope_sqtl < 0, pval_eqtl]
) %>%
  multiqq(flipY = TRUE)

# randomly selected genome-wide snps
qqplot_gw <- rbind(
  multiqq(list(GenomeWide = nomDF_genome[slope < 0, pval]))$data,
  multiqq(list(GenomeWide = nomDF_genome[slope > 0, pval]), flipY = T)$data
)

Title <- glue("{tissue}")
COLORS <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,6)], "grey")
names(COLORS) <- c("Productive", "Unproductive", "GenomeWide")


qqplot <- rbind(qqplot$data, qqplotNeg$data, qqplot_gw) %>%
  ggplot(aes(x, y, col = group)) + geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = Title) +
    scale_color_manual(values = COLORS) +
    theme_cowplot() +
    theme(legend.title = element_blank())


# save file

if (!dir.exists(out_prefix)) {
  dir.create(out_prefix, recursive = TRUE)
}

print("save qqplot")
ggsave(glue("{out_prefix}/{tissue}-qqplot.png"), qqplot, width = 7, height = 5, dpi = 200)

#---- plot enrichment ----

corr <- plotDF[, .(slope_sqtl, slope_eqtl, ctype = if_else(ctype == "PR", "PR", "UP"))]  %>% 
  split(by = "ctype") %>% 
  map(~cor.test(x = .x$slope_eqtl, y = .x$slope_sqtl, method = "p"))

corr.pvals <- map(corr, ~.x$p.value) %>% unlist
corr.estimates <- map(corr, ~.x$estimate[[1]]) %>% unlist

if (all(names(corr.pvals) == names(corr.estimates))) {
  corr.df <- data.frame(
    ctype = if_else(names(corr.pvals) == "PR", "Productive", "Unproductive"),
    pval = corr.pvals,
    estimate = corr.estimates
  )
  corr.df <- corr.df %>% 
  mutate(xpos = c(0, 0), ypos = c(0, 0))
} else {
  stop("names of pvals and estimates do not match")
}


scatter <- plotDF[, .(slope_sqtl, slope_eqtl, ctype)] %>% 
  mutate(ctype = if_else(ctype == "PR", "Productive", "Unproductive")) %>%
  ggplot() + geom_pointdensity(aes(slope_eqtl, slope_sqtl), alpha = .6) +
    scale_color_viridis_c() +
    geom_abline(aes(intercept = 0, slope = estimate), 
                data = corr.df, linetype = "dashed", color = "navy", linewidth = 1) +
    geom_text(
              aes(
                  x = xpos, y = ypos,
                  label = glue(
                               "cor: {corr}\np: {pvalue}", 
                               corr = if_else(abs(estimate) > .001, scales::number(estimate, .01), scales::scientific(estimate)),
                               pvalue = if_else(pval > .001, scales::number(pval, .01), scales::scientific(pval))
                              )
                  ),
              data = corr.df, size = 6, hjust = .5, vjust = 1
              ) +
    labs(x = "eQTL effect size", y = "sQTL effect size", title = glue("{tissue}")) +
    facet_wrap(~ctype) + 
    theme_cowplot() + 
    theme(strip.background = element_rect(fill = "white"))

print("save scatter")
ggsave(glue("{out_prefix}/{tissue}-scatter.png"), scatter, width = 7, height = 5, dpi = 200)

#---- save data ----
outobj <- list(
  joinedDF = mergeDF,
  corrDF = corr.df
)
print(glue("save data to {out_prefix}/{tissue}.rds"))
saveRDS(outobj, glue("{out_prefix}/{tissue}.rds"))

