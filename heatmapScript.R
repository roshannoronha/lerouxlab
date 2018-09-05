library(xlsx)
library(ggplot2)
library(pheatmap)
library(superheat)
library(RColorBrewer)
library(biomaRt)
biocLite("biomaRt")
setwd("~/LerouxLab/heatmap")

#read in heatmap datasheet
heatmap_df = read.csv("UpdatedHeatmapData.txt", sep="\t", stringsAsFactors = FALSE)
row_names = as.vector(heatmap_df[c(1:19848), 1])
heatmap_df = heatmap_df[c(1:19848), c(2:139)]
heatmap_df = sapply(heatmap_df, as.numeric)
rownames(heatmap_df) = row_names

allGenes_mat = pheatmap(heatmap_df, color = c("red", "yellow"), breaks =c(0, 0.5, 1), fontsize = 5, cellheight = 5, cellwidth = 5, border_color = "lightgrey", clustering_distance_rows = 'binary',cluster_cols = FALSE, clustering_distance_columns = 'binary', legend = TRUE, filename = "AllGenesHeatmap_Alphabetical.pdf")
write.csv(heatmap_df[allGenes_mat$tree_row$order,], quote = FALSE, file = "allGenesClustered2.csv")

#condensed heatmap
con_heatmap = read.csv("UpdatedHeatmapData.txt", sep="\t", stringsAsFactors = FALSE)
row_names = as.vector(con_heatmap[, 1])
con_heatmap = con_heatmap[, c(145:168)]
con_heatmap = apply(con_heatmap, 2, as.numeric)
rownames(con_heatmap) = row_names

#use this for now
con_mat = pheatmap(con_heatmap, clustering_distance_rows = 'binary',cluster_cols = FALSE, clustering_distance_columns = 'binary', color = colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))(10), treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, fontsize = 5, cellheight = 5, cellwidth = 5, border_color = "lightgrey", legend = TRUE, breaks = c(0,1,2, 5, 10, 20, 30, 40, 50, 60), legend_breaks = c(0,1,2, 5, 10, 20, 30, 40, 50, 60), file ="condensedHeatmapAllGenes_Alphabetical.pdf")
write.xlsx(con_heatmap[con_mat$tree_row$order,], "condensedgenesClustered.xlsx")


#adding leroux genes to allGenesHeatmapClustered
gene_df = read.csv("AllGenesHeatmapClusteredAnnotationMetadataOrthologues.csv", stringsAsFactors = FALSE)

lerouxGenes = read.csv("nrm.2017.60-s1.csv", stringsAsFactors = FALSE)
lerouxGenes = lerouxGenes[c(1:428),c(1:8)]


for (i in c(1:428)) {
  gene_df[which(lerouxGenes[i, 1] == gene_df[, 1]), c(2:10)] = c("CILIATED", paste(toString(lerouxGenes[i, 1]),"_", "CILIATED", sep = ""), lerouxGenes[i,c(2:8)])
}

write.table(gene_df, file = "AllGenesHeatmapUpdated_TEMP.csv", quote = FALSE, sep = "\t", col.names = NA)

#adding orthologues
#since some orthologues have multiple rows they will be saved in a seperate datasheet
#the order will match the clustered genes in AllGenesHeatmapClusteredAnnotationMetadataOrthologues

#read in list of genes from AllGenesHeatmapClusteredAnnotationMetadataOrthologues
allGenes = read.csv("Phylomaps/AllGenesHeatmapClusteredAnnotationMetadataORTHOLOGUES.csv", stringsAsFactors = FALSE)
geneNames = as.vector(allGenes[ c(6:nrow(allGenes)) , 1])

#read in orthologue datasheet
ortho = read.csv("orothologueData/celegansOrothologues.csv", stringsAsFactors = FALSE)

#empty dataframe to store orthologue data
storeClusteredOrtho = data.frame()

for (j in c(1:length(geneNames))) {
  storeClusteredOrtho = rbind(storeClusteredOrtho, ortho[which(geneNames[j] == ortho[, 5]), ])
}

write.table(storeClusteredOrtho, file = "ClusteredOrothologueData.csv", quote = FALSE, sep = "\t", col.names = NA)


###Get ensembl ID's#####
ensembl_hsapiens = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl_hsapiens)[1:5,]
attributes = listAttributes(ensembl_hsapiens)[1:5,]

#read in gene names from hsapiens_cildb
hSapiens_cilDB = read.csv("cilDBData/hsapiens_cildb.csv", stringsAsFactors = FALSE)
#read in gene names from phylomap
phylomap = read.csv("Phylomaps/AllGenesHeatmapClusteredAnnotationMetadataCILIATEDNOTCILIATEDSORTED.csv", stringsAsFactors = FALSE)
phylomap = phylomap[c(6:nrow(phylomap)), ]
phylomapGeneName = phylomap[,1]

#get ensembl id's based on gene names
ensemblIds = getBM(attributes= c("external_gene_name", "ensembl_gene_id"), filter = "external_gene_name", values = phylomapGeneName, mart = ensembl_hsapiens)

#get sum and store in the appropriate row and column
for (i in c(1:nrow(ensemblIds))) {
  sums = c(apply(hSapiens_cilDB[which(hSapiens_cilDB[,4] == ensemblIds[i,2]),][, c(10,11,12)], 2, sum), ensemblIds[i,2])
  phylomap[which(phylomap[, 1] == ensemblIds[i,1]), c(3,4,5, 14)] = sums
}

write.table(phylomap, file = "AllGenesHeatmapClusteredAnnotationMetadataSorted_EVIDENCE.csv", quote = FALSE, sep = "\t")

