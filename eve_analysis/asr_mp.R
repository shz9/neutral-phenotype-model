library(ape)
library(phangorn)
library(rjson)


#args <- commandArgs(trailingOnly=TRUE)
args <- c("~/10:100703223-100704293.json")
output_dir <- "~/"

json_data <- fromJSON(file = args[1])

tree <- read.tree(text=json_data$tree)

# Clean up the tip names in the tree:
for (tl in 1:length(tree$tip.label)){
  tree$tip.label[tl] <- sub('(_[^_]+){3}$', '', tree$tip.label[tl])
  tree$tip.label[tl] <- sub('_scaffold', '', tree$tip.label[tl])
}

seqs <- do.call("rbind", lapply(json_data$alignments, as.data.frame))
seqs$species <- gsub("\\s*\\[[^\\)]+\\]", "", seqs$species)

tip_seqs <- seqs[seqs$species %in% tree$tip.label, c('species', 'seq')]

con <- textConnection(paste(paste0('>', tip_seqs$species), 
                            tip_seqs$seq, sep="\n"))
phydat <- read.phyDat(con, format="fasta", type="USER",
                      levels=c("a","c","g","t","-"),
                      ambiguity = c("n"))
close(con)

anc_pars <- ancestral.pars(tree, phydat, return="phyDat")

for (i in 1:length(anc_pars)){
  
  if (!names(anc_pars)[i] %in% tree$tip.label){
    names(anc_pars)[i] <- tree$node.label[as.integer(names(anc_pars)[i])-Ntip(tree)]
  }
  
}

#char_anc <- as.character(anc_pars)

#for (n in names(anc_pars)){
#  pars_seq <- char_anc[n,]
#  ml_seq <- strsplit(seqs[seqs$species == n, 'seq'], "")[[1]]
#  
#  print(n)
#  print(sum(toupper(pars_seq) == ml_seq) / length(pars_seq)) 
#}


write.phyDat(anc_pars, paste0(output_dir, 
                              gsub("\\.json", ".fa", basename(args[1]))
                              ), format = "fasta")



