
library(GOplot)
library(topGO)
library(ggplot2)
setwd("/vol/storage/data/GT2021/GOPM/")
universe <- readMappings("universe") 
universe_genelist <- names(universe)
filename_tuq1 = "CDScpgtopgo"
filename_tuq2 = "CDSchhtopgo"
filename_tuq3 = "CDSchgtopgo"
testset_genelist <- readLines(filename_tuq1)
gene_list <- factor(as.integer(universe_genelist %in% testset_genelist))
names(gene_list) <- universe_genelist
ontology_type = "BP"
GO_data <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list, annot = annFUN.gene2GO, gene2GO=universe)

result_topGO_w01 <- runTest(GO_data, algorithm = "weight01", statistic = "fisher")
result_table_w01 <- GenTable(GO_data, Fisher = result_topGO_w01, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 

result_topGO_pc <- runTest(GO_data, algorithm = "parentchild", statistic = "fisher")
result_table_pc <- GenTable(GO_data, Fisher = result_topGO_pc, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 
pdf(file = paste0("topGO_result_PC_significant_", filename_tuq1, "_", ontology_type, ".pdf"))
par(cex=0.2)
showSigOfNodes(GO_data, score(result_topGO_pc), firstSigNodes = 10, useInfo = "def")
dev.off()
pdf(file = paste0("topGO_result_weight01_significant_", filename_tuq1, "_", ontology_type, ".pdf"))
showSigOfNodes(GO_data, score(result_topGO_w01), firstSigNodes = 10, useInfo = "def")
dev.off()
short_res <-as.data.frame(result_table_pc[,c(2,4)])
res_short<-result_table_pc[which(as.integer(result_table_pc$Fisher)<0.05),]
short_res <-as.data.frame(res_short[,c(2,4)])
pdf(file = paste0("topGO_result_PC_hist_significant_", filename_tuq1, "_", ontology_type, ".pdf"))
par(mfrow = c(1, 2))
ggplot(data=short_res, aes(x=Term, y=Significant, fill=Term)) +
  geom_bar(width=0.8, stat = "identity") +
  theme(axis.text.x=element_text(face = "bold", lineheight = 20, size = 10,angle = 70, hjust = 1))
dev.off()

write.csv(result_table_pc, file=paste0("topGO_result_parentchild_", filename_tuq1, "_", ontology_type, ".csv"))
write.table(result_table_pc[which( as.integer(result_table_pc$Fisher)<0.05),], file=paste0("topGO_result_parentchild_significant_", filename_tuq1, "_", ontology_type,".csv"))
write.csv(result_table_w01, file=paste0("topGO_result_weight01_", filename_tuq1, "_", ontology_type, ".csv"))
write.csv(result_table_w01[which(result_table_w01$Fisher<0.05),], file=paste0("topGO_result_weight01_significant_", filename_tuq1, "_", ontology_type, ".csv"))



testset_genelist <- readLines(filename_tuq2)
gene_list <- factor(as.integer(universe_genelist %in% testset_genelist))
names(gene_list) <- universe_genelist
ontology_type = "BP"
GO_data <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list, annot = annFUN.gene2GO, gene2GO=universe)

result_topGO_w01 <- runTest(GO_data, algorithm = "weight01", statistic = "fisher")
result_table_w01 <- GenTable(GO_data, Fisher = result_topGO_w01, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 

result_topGO_pc <- runTest(GO_data, algorithm = "parentchild", statistic = "fisher")
result_table_pc <- GenTable(GO_data, Fisher = result_topGO_pc, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 
pdf(file = paste0("topGO_result_PC_significant_", filename_tuq2, "_", ontology_type, ".pdf"))
par(cex=0.2)
showSigOfNodes(GO_data, score(result_topGO_pc), firstSigNodes = 10, useInfo = "def")
dev.off()
pdf(file = paste0("topGO_result_weight01_significant_", filename_tuq2, "_", ontology_type, ".pdf"))
showSigOfNodes(GO_data, score(result_topGO_w01), firstSigNodes = 10, useInfo = "def")
dev.off()

short_res <-as.data.frame(result_table_pc[,c(2,4)])
res_short<-result_table_pc[which(as.integer(result_table_pc$Fisher)<0.05),]
short_res <-as.data.frame(res_short[,c(2,4)])
pdf(file = paste0("topGO_result_PC_hist_significant_", filename_tuq2, "_", ontology_type, ".pdf"))
par(mfrow = c(1, 2))
ggplot(data=short_res, aes(x=Term, y=Significant, fill=Term)) +
  geom_bar(width=0.8, stat = "identity") +
  theme(axis.text.x=element_text(face = "bold", lineheight = 20, size = 10,angle = 70, hjust = 1))
dev.off()
write.csv(result_table_pc, file=paste0("topGO_result_parentchild_", filename_tuq2, "_", ontology_type, ".csv"))
write.table(result_table_pc[which( as.integer(result_table_pc$Fisher)<0.05),], file=paste0("topGO_result_parentchild_significant_", filename_tuq2, "_", ontology_type,".csv"))
write.csv(result_table_w01, file=paste0("topGO_result_weight01_", filename_tuq2, "_", ontology_type, ".csv"))
write.csv(result_table_w01[which(result_table_w01$Fisher<0.05),], file=paste0("topGO_result_weight01_significant_", filename_tuq2, "_", ontology_type, ".csv"))

testset_genelist <- readLines(filename_tuq3)
gene_list <- factor(as.integer(universe_genelist %in% testset_genelist))
names(gene_list) <- universe_genelist
ontology_type = "BP"
GO_data <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list, annot = annFUN.gene2GO, gene2GO=universe)

result_topGO_w01 <- runTest(GO_data, algorithm = "weight01", statistic = "fisher")
result_table_w01 <- GenTable(GO_data, Fisher = result_topGO_w01, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 

result_topGO_pc <- runTest(GO_data, algorithm = "parentchild", statistic = "fisher")
result_table_pc <- GenTable(GO_data, Fisher = result_topGO_pc, orderBy = "Fisher", ranksOf = "Fisher", topNodes = 1000) 
pdf(file = paste0("topGO_result_PC_significant_", filename_tuq3, "_", ontology_type, ".pdf"))
par(cex=0.2)
showSigOfNodes(GO_data, score(result_topGO_pc), firstSigNodes = 10, useInfo = "def")
dev.off()
pdf(file = paste0("topGO_result_weight01_significant_", filename_tuq3, "_", ontology_type, ".pdf"))
showSigOfNodes(GO_data, score(result_topGO_w01), firstSigNodes = 10, useInfo = "def")
dev.off()

short_res <-as.data.frame(result_table_pc[,c(2,4)])
res_short<-result_table_pc[which(as.integer(result_table_pc$Fisher)<0.05),]
short_res <-as.data.frame(res_short[,c(2,4)])
pdf(file = paste0("topGO_result_PC_hist_significant_", filename_tuq3, "_", ontology_type, ".pdf"))
par(mfrow = c(1, 2))
ggplot(data=short_res, aes(x=Term, y=Significant, fill=Term)) +
  geom_bar(width=0.8, stat = "identity") +
  theme(axis.text.x=element_text(face = "bold", lineheight = 20, size = 10,angle = 70, hjust = 1))
dev.off()

write.csv(result_table_pc, file=paste0("topGO_result_parentchild_", filename_tuq3, "_", ontology_type, ".csv"))
write.csv(result_table_pc[which( as.integer(result_table_pc$Fisher)<0.05),], file=paste0("topGO_result_parentchild_significant_", filename_tuq3, "_", ontology_type,".csv"))
write.csv(result_table_w01, file=paste0("topGO_result_weight01_", filename_tuq3, "_", ontology_type, ".csv"))
write.csv(result_table_w01[which(result_table_w01$Fisher<0.05),], file=paste0("topGO_result_weight01_significant_", filename_tuq3, "_", ontology_type, ".csv"))
