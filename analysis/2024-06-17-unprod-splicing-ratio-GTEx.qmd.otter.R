












#| label: setup
#| include: false
#| cache: false
library(tidyverse)
library(data.table)
library(glue)


library(cowplot)
theme_set(theme_cowplot())


library(furrr)
plan(multisession, workers = 6)






tissues <- dir("/project/yangili1/cdai/SpliFi/code/results/pheno/noisy/GTEx")
tissues <- tissues[!str_detect(tissues, "Bladder")]





computeUnprodRatio <- function(tissue_rds) {

    introns <- readRDS(tissue_rds)
    datacols <- names(introns)[5:ncol(introns)]
    # remove clu_type == N, clusters with only unproductive introns
    introns  <- introns[str_detect(clu_type, 'PR')]
    # total reads (by sample)
    totals <- colSums(introns[, ..datacols])
    # unproudctive reads (by sample)
    unprod <- introns[intron_type == 'UP', ..datacols] %>% colSums
    # unproductive reads / total reads
    if (all(names(totals) == names(unprod))) {
        unprod_ratio <- round(unprod / totals, 5)
    } else {
        stop('samples of totals and unprod are not the same')
    }

    return(unprod_ratio)

}


point_scale <- function(x, newmin, newmax) {
  newmax  <- max(newmax)
  newmin  <- min(newmin)
  new_range <- abs(newmax - newmin) 
  old_range <- abs(max(x) - min(x))

  # Normalize the values between 0 and 1
  normalized_x <- (x - min(x)) / old_range

  # Scale the normalized values to the new range and offset by 0.5
  scaled_x <- normalized_x * new_range + newmin

  return(scaled_x)
}






#| label: read-splice_ratio_by_tissue_rds

tissues_rds <- glue("/project/yangili1/cdai/SpliFi/data/ExtractFractions/GTEx/{tissues}.numerators_constcounts.noise_by_intron.rds")
names(tissues_rds) <- tissues




#| label: "read all gtex ratio rds files"
ratios <- map(tissues_rds, computeUnprodRatio)




ratios <- imap_dfr(
    ratios,
    \(x, y) enframe(x, "sample", "unprodRatio") %>% add_column(tissue = y)
)
ratios  <- as.data.table(ratios) %>%
    .[, .(sample, unprodRatio, nsample = uniqueN(sample)), by = tissue]





ratios[, .(median_prod_ratio = median(unprodRatio)), by = tissue][order(-median_prod_ratio)]

median_unprod_ratio_all_gtex <- ratios[tissue != "EUR", .(median_prod_ratio = median(unprodRatio)), by = tissue
    ][order(-median_prod_ratio)] %>% 
    pull(median_prod_ratio) %>% median




tissue.labels <- fread("/project/yangili1/cdai/splice-pub/analysis/gtex-tissue-code.csv",
                       header = F, col.names = c("tissue", "label"))




df <- ratios[, .(sample, unprodRatio, nsample, med = median(unprodRatio),rk=rank(unprodRatio, ties.method = "first")), by = tissue
        ][, .(sample, unprodRatio, nsample, med, tissue = forcats::fct_reorder(tissue, -med), rk)
        ][, .(sample, tissue, unprodRatio, nsample, med, tissuei = as.integer(tissue), rk)
        ][, .(sample, unprodRatio, nsample, med, tissuei, rk, pointx = point_scale(rk, (tissuei - .5), tissuei+.5)), by = tissue
        ][order(tissuei, rk)] %>% 
    left_join(y = tissue.labels, by = "tissue")






labels  <- df[, .(tissue, label, tissuei, nsample)] %>% unique
labels




#| label: fig-unprod_ratio_by_tissue
#| fig-cap: "Unproductive splicing ratio in GTEx across tissues"
#| fig-width: 10
#| fig-height: 5
#| out-width: 100%

label_x1 <- function(ids) {
    labels[tissuei == ids, label]
}

label_x2 <- function(ids) {
    labels[tissuei == ids, nsample] 
}


df %>% left_join(y = tissue.labels, by = "tissue") %>%
      ggplot() +
        geom_point(aes(x = pointx, y = unprodRatio), color = 'midnightblue', alpha = .1) +
        geom_tile(aes(x = tissuei, y = med), color = 'firebrick', fill= 'firebrick', height = 5e-5, width = .9) +
        scale_x_continuous(
            breaks = unique(df$tissuei), 
            labels = unique(df$label),
            expand = c(0, .5),
            guide = guide_axis(angle = 90),
            sec.axis = sec_axis(~., 
                                breaks = unique(df$tissuei),
                                labels = label_x2,  
                                guide = guide_axis(angle = 90)
                               )
            ) +
        scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
        labs(x = "Tissue / Number of Samples",
            y = "Fraction of unproductive reads"
            )






#| include: false
#| eval: false
#zzzz
median_unprod_ratio_all_gtex 

df[str_detect(tissue, "Testis")] %>% head





#| include: false
#| eval: false

httpgd::hgd(host = "10.50.250.37", port = 9002, token = F)

