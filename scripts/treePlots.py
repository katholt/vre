# with mlst
t<-read.tree("~/code/vre/Phylogenetics/final_trees/95_gene_cons_RAxML_bipartitions.rooted.tree")
mlst<-"~/code/vre/All_MLST_Results.csv"
source("~/code/holtlab/Rcode/plotTree.R")
plotTree(t,infoFile=mlst,tip.labels=T)

