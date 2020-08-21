library(ape)
library(data.table)
library(ggplot2)
library(geiger)


generate_aicc_bootstrap <- function(sp_tree, trait_tree_df, trait_name) {

    tdf <- as.data.frame(lapply(trait_tree_df[trait_name], as.numeric))
    rownames(tdf) <- rownames(trait_tree_df)
    bm_res <- fitContinuous(sp_tree, tdf, model="BM")
    ou_res <- fitContinuous(sp_tree, tdf, model="OU")
    return(list(trait_name,
                as.numeric(bm_res$opt$aicc),
                as.numeric(ou_res$opt$aicc)))
}

main_input_dirs <- list("Real Sequences" = "./data/2_q_traits/0_chopped_alignments",
                        "Simulated Sequences" = "./data/2_q_traits/2_simulated_sequences")
tree_dir <- "./metadata/segment_trees/paml_inferred/"

plot_trait_names <- list("A_percent" = "Percent of A",
                         "gc_content" = "GC Content",
                         "longest_ORF" = "Longest ORF")

output_dir <- "./plots/geiger_inference/"


for (dir_type in names(main_input_dirs)) {
    
    print(dir_type)
    parent_dir <- main_input_dirs[[dir_type]]
    
    final_df <- data.frame(trait=character(),
                           BM_mean_aicc=numeric(),
                           OU_mean_aicc=numeric(),
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
            final_df[nrow(final_df) + 1, ] <- generate_aicc_bootstrap(species_tree,
                                                                      t_trait_df,
                                                                      tr)
            final_df[nrow(final_df), "trait"] <- plot_trait_names[[final_df[nrow(final_df), "trait"]]]
        }
    }

    if (nrow(final_df) > 0) {
        
        final_df$trait <- as.factor(final_df$trait)

        aic_f_df <- final_df[, 2:ncol(final_df)]

        min_aicc <- min(aic_f_df) - .05 * sign(min(aic_f_df)) * min(aic_f_df)
        max_aicc <- max(aic_f_df) + .05 * sign(max(aic_f_df)) * max(aic_f_df)

        top_polygon_df <- data.frame(x = c(min_aicc, min_aicc, max_aicc),
                                     y = c(min_aicc, max_aicc, max_aicc))

        bottom_polygon_df <- data.frame(x=c(min_aicc, max_aicc, max_aicc),
                                        y=c(min_aicc, min_aicc, max_aicc))

        p <- ggplot(final_df, aes(x = BM_mean_aicc, y = OU_mean_aicc)) 
        p <- p + geom_point(aes(color = trait))
        p <- p + geom_polygon(data = top_polygon_df, aes(x=x, y=y), fill="#F4CC7060")
        p <- p + geom_polygon(data = bottom_polygon_df, aes(x=x, y=y), fill="#6AB18760")
        p <- p + labs(x = "BM corrected AIC",
                      y = "OU corrected AIC",
                      title = dir_type)
        p <- p + annotate("text", x = min(aic_f_df), y = max(aic_f_df), label = "BM favored",
                          vjust="inward", hjust="inward")
        p <- p + annotate("text", x = max(aic_f_df), y = min(aic_f_df), label = "OU favored",
                          vjust="inward", hjust="inward")

        ggsave(paste0(output_dir, dir_type, ".pdf"))
    }
}

