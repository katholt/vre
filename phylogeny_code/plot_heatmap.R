source('vlsci/VR0082/shared/code/holtlab/getRecomb.R')
source('/vlsci/VR0082/shared/code/holtlab/Rcode/plotTree.R')

tree.file <- read.tree('/vlsci/VR0082/shared/data/enterococcus/faecium/RedDogv51_How2013Alf2015_CP006620/CP006620_CP006620_alleles_var_cons0.95.tree')
ingroup.tree<-drop.tip(tree.file,c('aus0013_AC026VACXX_CAGATC_L005','aus0028_AC026VACXX_ATCACG_L006','aus0103_AC026VACXX_GATCAG_L008'))
recombination.density <- getRecombHeatmap(snp.table,ref="Ref",w=10000)
plotTree(tree=ingroup.tree, heatmapData=recombination.density,heatmap.colours=colorRampPalette(c("white","orange","red"),space="rgb")(100), colLabelCex=0.0001)
