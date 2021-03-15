library(ouch)
library(ape)
library(data.table)
library(ggplot2)

########### Global variables ###########

bootstrap_res <- F
main_input_dirs <- list("Real Sequences" = "./data/3_q_traits/0_chopped_alignments",
                        "Simulated Sequences" = "./data/3_q_traits/2_simulated_sequences")
tree_dir <- "./metadata/segment_trees/paml_inferred/"

plot_trait_names <- list("A_percent" = "Percent of A",
                         "gc_content" = "GC Content",
                         "longest_ORF" = "Longest ORF",
                         "longest_ORF_alan" = "Longest ORF (Alan)")

output_dir <- "./ouch_inference/"

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

generate_bm_vs_ou_results <- function(sp_ouch_tree, trait_tree_df, segment, trait_name) {

    bm_res <- brown(as.data.frame(lapply(trait_tree_df[trait_name], as.numeric)),
                    sp_ouch_tree)
    bm_res <- summary(bm_res)

    hh_res <- hansen(tree=sp_ouch_tree,
                     data=as.data.frame(lapply(trait_tree_df[trait_name], as.numeric)),
                     regimes=trait_tree_df["regimes"],
                     sqrt.alpha=1,
                     sigma=1,
                     maxit=1000,
                     fit=T)
    hh_res <- summary(hh_res)

    return(list(segment,
                trait_name,
                bm_res$loglik,
                bm_res$aic.c,
                unlist(bm_res$theta),
                hh_res$loglik,
                hh_res$aic.c,
                unlist(hh_res$optima)))
}

########### Running the analysis ###########

setwd(dirname(parent.frame(2)$ofile))
dir.create(output_dir, recursive = T)

for (dir_type in names(main_input_dirs)) {
    
    print(dir_type)
    parent_dir <- main_input_dirs[[dir_type]]

    if (bootstrap_res) {
        final_df <- data.frame(segment=character(),
                           trait=character(),
                           BM_mean_aicc=numeric(),
                           BM_hi_ci_aicc=numeric(),
                           BM_lo_ci_aicc=numeric(),
                           OU_mean_aicc=numeric(),
                           OU_hi_ci_aicc=numeric(),
                           OU_lo_ci_aicc=numeric(),
                           stringsAsFactors=F)
    }
    else{
        final_df <- data.frame(segment=character(),
                           trait=character(),
                           BM_loglik=numeric(),
                           BM_aic_c=numeric(),
                           BM_Zeq=numeric(),
                           OU_loglik=numeric(),
                           OU_aic_c=numeric(),
                           OU_Zeq=numeric(),
                           stringsAsFactors=F)
    }

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
            if (bootstrap_res) {
                final_df[nrow(final_df) + 1, ] <- generate_aicc_bootstrap(sp_ouch_tree,
                                                                          trait_tree_df,
                                                                          segment_name,
                                                                          tr)
            }
            else{
                final_df[nrow(final_df) + 1, ] <- generate_bm_vs_ou_results(sp_ouch_tree,
                                                                          trait_tree_df,
                                                                          segment_name,
                                                                          tr)
            }
            final_df[nrow(final_df), "trait"] <- plot_trait_names[[final_df[nrow(final_df), "trait"]]]
        }
    }

    if (nrow(final_df) > 0) {
        write.csv(final_df, paste0(output_dir, dir_type, ".csv"))
    }
}

