library(ape)
library(data.table)
library(ggplot2)
library(geiger)


generate_bm_vs_ou_results_geiger <- function(sp_tree, trait_tree_df, segment, trait_name) {

    tdf <- as.data.frame(lapply(trait_tree_df[trait_name], as.numeric))
    rownames(tdf) <- rownames(trait_tree_df)
    
    bm_res <- fitContinuous(sp_tree, tdf, model="BM",
                            control=list(method=c("L-BFGS-B"), 
                                         niter=100, FAIL=1e200, 
                                         hessian=FALSE, 
                                         CI=0.95),
                            bounds=list(sigsq=c(min=1e-12,max=1e12))
    )
    ou_res <- fitContinuous(sp_tree, tdf, model="OU",
                            control=list(method=c("L-BFGS-B"),
                                         niter=100, FAIL=1e200, 
                                         hessian=FALSE, CI=0.95),
                            bounds=list(sigsq=c(min=1e-12,max=1e12), 
                                        alpha=c(min=1e-12,max=1e12))
                            )
    
    return(list(segment,
                trait_name,
                as.numeric(bm_res$opt$lnL),
                as.numeric(bm_res$opt$aicc),
                as.numeric(bm_res$opt$z0),
                as.numeric(ou_res$opt$lnL),
                as.numeric(ou_res$opt$aicc),
                as.numeric(ou_res$opt$z0)
                )
           )
}

main_input_dirs <- list("Real Sequences" = "./data/3_q_traits/0_chopped_alignments",
                        "Simulated Sequences" = "./data/3_q_traits/2_simulated_sequences")
tree_dir <- "./metadata/segment_trees/paml_inferred/"

plot_trait_names <- list("A_percent" = "Percent of A",
                         "gc_content" = "GC Content",
                         "longest_ORF" = "Longest ORF",
                         "longest_ORF_alan" = "Longest ORF (Alan)")

output_dir <- "./geiger_inference/"

########### Running the analysis ###########

setwd(dirname(parent.frame(2)$ofile))
dir.create(output_dir, recursive = T)

for (dir_type in names(main_input_dirs)) {
    
    print(dir_type)
    parent_dir <- main_input_dirs[[dir_type]]
    
    final_df <- data.frame(segment=character(),
                           trait=character(),
                           BM_loglik=numeric(),
                           BM_aicc=numeric(),
                           BM_z0=numeric(),
                           OU_loglik=numeric(),
                           OU_aicc=numeric(),
                           OU_z0=numeric(),
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
        rownames(t_trait_df) <- colnames(trait_df)

        for (tr in rownames(trait_df)) {
            print(tr)
            final_df[nrow(final_df) + 1, ] <- generate_bm_vs_ou_results_geiger(species_tree,
                                                                      t_trait_df,
                                                                      segment_name,
                                                                      tr)
            final_df[nrow(final_df), "trait"] <- plot_trait_names[[final_df[nrow(final_df), "trait"]]]
        }
    }

    if (nrow(final_df) > 0) {
        write.csv(final_df, paste0(output_dir, dir_type, ".csv"))
    }
}

