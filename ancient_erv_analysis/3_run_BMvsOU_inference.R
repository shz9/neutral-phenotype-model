library(ouch)
library(ape)
library(data.table)
library(ggplot2)

########### Global variables ###########

main_input_dirs <- list("Real Sequences" = "./data/2_q_traits/0_chopped_alignments",
                        "Simulated Sequences" = "./data/2_q_traits/2_simulated_sequences")
tree_dir <- "./metadata/segment_trees/paml_inferred/"

plot_trait_names <- list("A_percent" = "Percent of A",
                         "gc_content" = "GC Content",
                         "longest_ORF" = "Longest ORF",
                         "longest_ORF_alan" = "Longest ORF (Alan)")

plot_dir <- "./plots/ouch_inference/"
output_dir <- "./tables/ouch_inference/"

########### Helper functions ###########

generate_aicc_bootstrap <- function(sp_ouch_tree, trait_tree_df, segment, trait_name) {
    
    bm_res <- brown(as.data.frame(lapply(trait_tree_df[trait_name], as.numeric)),
                    sp_ouch_tree)
    bm_boot_res <- bootstrap(bm_res)
    
    hh_res <- hansen(tree=sp_ouch_tree,
                     data=as.data.frame(lapply(trait_tree_df[trait_name], as.numeric)),
                     regimes=trait_tree_df["regimes"],
                     sqrt.alpha=1,
                     sigma=1,
                     maxit=1000,
                     fit=T)
    hh_boot_res <- bootstrap(hh_res)
    return(list(segment,
                trait_name,
                mean(bm_boot_res$aic.c),
                mean(bm_boot_res$aic.c) + 1.96*sd(bm_boot_res$aic.c),
                mean(bm_boot_res$aic.c) - 1.96*sd(bm_boot_res$aic.c),
                mean(hh_boot_res$aic.c),
                mean(hh_boot_res$aic.c) + 1.96*sd(hh_boot_res$aic.c),
                mean(hh_boot_res$aic.c) - 1.96*sd(hh_boot_res$aic.c)))
}

########### Running the analysis ###########

setwd(dirname(parent.frame(2)$ofile))
dir.create(plot_dir, recursive = T)
dir.create(output_dir, recursive = T)

for (dir_type in names(main_input_dirs)) {
    
    print(dir_type)
    parent_dir <- main_input_dirs[[dir_type]]
    
    final_df <- data.frame(segment=character(),
                           trait=character(),
                           BM_mean_aicc=numeric(),
                           BM_hi_ci_aicc=numeric(),
                           BM_lo_ci_aicc=numeric(),
                           OU_mean_aicc=numeric(),
                           OU_hi_ci_aicc=numeric(),
                           OU_lo_ci_aicc=numeric(),
                           stringsAsFactors=F)

    for (input_dir in list.dirs(parent_dir, recursive=F)) {
        
        print(input_dir)
        segment_name <- basename(input_dir)
        # Read the tree and trait data
        species_tree <- read.tree(paste0(tree_dir, segment_name, ".nwk"))
        trait_df <- read.csv(paste0(input_dir, "/traits.csv"), row.names=1)
        
        for (t in species_tree$tip.label) {
            if (!(t %in% colnames(trait_df))) {
                species_tree <- drop.tip(species_tree, t)
            }
        }

        # Transform the trait data
        trait_df <- trait_df[, c(species_tree$tip.label)]
        t_trait_df <- transpose(trait_df)
        colnames(t_trait_df) <- rownames(trait_df)
        t_trait_df$species <- colnames(trait_df)

        sp_ouch_tree <- ape2ouch(species_tree)
        sp_ot_df <- as(sp_ouch_tree, "data.frame")

        trait_tree_df <- merge(sp_ot_df, t_trait_df,
                               by.x="labels", by.y="species", all.x=T)
        trait_tree_df$regimes <- as.factor("global")

        for (tr in rownames(trait_df)) {
            print(tr)
            final_df[nrow(final_df) + 1, ] <- generate_aicc_bootstrap(sp_ouch_tree,
                                                                      trait_tree_df,
                                                                      segment_name,
                                                                      tr)
            final_df[nrow(final_df), "trait"] <- plot_trait_names[[final_df[nrow(final_df), "trait"]]]
        }
    }

    if (nrow(final_df) > 0) {
        
        write.csv(final_df, paste0(output_dir, dir_type, ".csv"))
        final_df$trait <- as.factor(final_df$trait)
        final_df <- within(final_df, rm("segment"))

        aic_f_df <- final_df[, 2:ncol(final_df)]

        min_aicc <- min(aic_f_df) - .05 * sign(min(aic_f_df)) * min(aic_f_df)
        max_aicc <- max(aic_f_df) + .05 * sign(max(aic_f_df)) * max(aic_f_df)

        top_polygon_df <- data.frame(x = c(min_aicc, min_aicc, max_aicc),
                                     y = c(min_aicc, max_aicc, max_aicc))

        bottom_polygon_df <- data.frame(x=c(min_aicc, max_aicc, max_aicc),
                                        y=c(min_aicc, min_aicc, max_aicc))

        p <- ggplot(final_df, aes(x = BM_mean_aicc, y = OU_mean_aicc)) 
        p <- p + geom_point(aes(color = trait))
        p <- p + geom_errorbar(aes(ymin = OU_lo_ci_aicc, ymax = OU_hi_ci_aicc, color = trait))
        p <- p + geom_errorbarh(aes(xmin = BM_lo_ci_aicc, xmax = BM_hi_ci_aicc, color = trait))
        p <- p + geom_polygon(data = top_polygon_df, aes(x=x, y=y), fill="#F4CC7060")
        p <- p + geom_polygon(data = bottom_polygon_df, aes(x=x, y=y), fill="#6AB18760")
        p <- p + labs(x = "BM corrected AIC",
                      y = "OU corrected AIC",
                      title = dir_type)
        p <- p + annotate("text", x = min(aic_f_df), y = max(aic_f_df), label = "BM favored",
                          vjust="inward", hjust="inward")
        p <- p + annotate("text", x = max(aic_f_df), y = min(aic_f_df), label = "OU favored",
                          vjust="inward", hjust="inward")

        ggsave(paste0(plot_dir, dir_type, ".pdf"))
    }
}

