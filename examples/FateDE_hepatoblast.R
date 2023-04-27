library(monocle)
library(FateDE)

data(hepatoblast)

#### Trajectory inference ####

fData <- data.frame(gene_short_name = hepatoblast$feature_info$feature_id, row.names = hepatoblast$feature_info$feature_id)
pData <- as.data.frame(hepatoblast$cell_info)
rownames(pData) = pData$cell_id

pd <- new("AnnotatedDataFrame", data = pData)
fd <- new("AnnotatedDataFrame", data = fData)
expr = t(hepatoblast$expression)

HSMM <- newCellDataSet(as.matrix(expr),
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))


######### Clustering cells without marker genes ###########
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 3,
                        reduction_method = 'tSNE', verbose = T)

HSMM <- clusterCells(HSMM, num_clusters = 6)
plot_cell_clusters(HSMM, 1, 2, color = "day")

######### Constructing Single Cell Trajectories #########
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~Cluster")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)

HSMM <- reduceDimension(HSMM,reduction_method = 'DDRTree')

HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "Cluster")
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "milestone_id")
plot_cell_trajectory(HSMM, color_by = "putative_cell_type")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "Cluster")+ geom_point(color = 'gray') + theme(legend.position = 'none')+
  theme(text = element_blank(), axis.ticks.length = unit(0,'cm'),axis.line.x=element_blank(),axis.line.y=element_blank())

######### Differential Expression Analysis #########

BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


p=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,pval < 1e-1)),],
                              branch_point = 1,
                              num_clusters = 3,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T,
                              return_heatmap=TRUE)

row_cluster <- cutree(p$ph_res$tree_row,k=3)
G1=names(row_cluster[row_cluster %in% c(2)])
G2=names(row_cluster[row_cluster %in% c(1)])

##### fitting #######
expr = data.frame(t(as.matrix(HSMM@assayData[["exprs"]]) ))
colnames(expr) = HSMM@featureData@data[["gene_short_name"]]

tmp = data.frame(
  State = HSMM@phenoData@data[["State"]],
  Pseudotime = HSMM@phenoData@data[["Pseudotime"]]
)

rownames(tmp) = tmp$cell_id
tmp = cbind(tmp,expr)


df1_notfit = subset(tmp, milestone_id != 3)
df1_notfit = df1_notfit[order(df1_notfit$Pseudotime),]
df2_notfit = subset(tmp, milestone_id != 2)
df2_notfit = df2_notfit[order(df2_notfit$Pseudotime),]

df1_fit = fateDE_fit_one(df1_notfit)
df2_fit = fateDE_fit_one(df2_notfit)

##### FateDE #######
result = FateDE_main(df1_fit, df2_fit, G1, G2, GP_FSV = 0.6)
G1_1 = result$G1_1
G2_1 = result$G2_1

#################### plot ###########################

# single gene curve
gene = 'Nr1h4'
plot_SingleGene_trend(list(df1_notfit,df2_notfit),df1_fit,df2_fit,gene)+
  scale_color_manual(values=c("#999999", "#0073C2FF",'#c74545'))

# double genes curve
double_gene = c('Zbtb8a','Zfp791')
plot_DoubleGenes_trend(list(df1_fit,df2_fit),
                       df1_fit,
                       df2_fit,
                       double_gene)

# heatmap
plot_Genes_Heatmap(df1_fit,df2_fit,G1_1)
plot_Genes_Heatmap(df1_fit,df2_fit,G2_1)

#################### Enrich Analysis ###########################
# 1.load
library(AnnotationDbi)
library(clusterProfiler)
library(KEGG.db)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

# 2.Gene ID conversion
gene.df = bitr(c(G1_1,G2_1), fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db)
gene = gene.df$ENTREZID

# 3.GO Enrich
ego_ALL = enrichGO(gene = gene,
                   OrgDb = org.Mm.eg.db,
                   keyType = 'ENTREZID',
                   ont = 'ALL',
                   pAdjustMethod = 'BH',
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_CC = enrichGO(gene = gene,
                  OrgDb = org.Mm.eg.db,
                  keyType = 'ENTREZID',
                  ont = 'CC',
                  pAdjustMethod = 'BH',
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)


ego_BP = enrichGO(gene = gene,
                  OrgDb = org.Mm.eg.db,
                  keyType = 'ENTREZID',
                  ont = 'BP',
                  pAdjustMethod = 'BH',
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

ego_MF = enrichGO(gene = gene,
                  OrgDb = org.Mm.eg.db,
                  keyType = 'ENTREZID',
                  ont = 'MF',
                  pAdjustMethod = 'BH',
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# 4.Filter some pathways and combine them
display_number = c(22,10,10)  # Number of selected paths
ego_result_BP = as.data.frame(ego_BP)
ego_result_BP = ego_result_BP[order(ego_result_BP$Count,decreasing = TRUE),][1:display_number[1],]

ego_result_CC = as.data.frame(ego_CC)
ego_result_CC = ego_result_CC[order(ego_result_CC$Count,decreasing = TRUE),][1:display_number[2],]

ego_result_MF = as.data.frame(ego_MF)
ego_result_MF = ego_result_MF[order(ego_result_MF$Count,decreasing = TRUE),][1:display_number[3],]


go_enrich_df = data.frame(
  ID = c(ego_result_BP$ID,ego_result_CC$ID,ego_result_MF$ID),
  GeneNumber = c(ego_result_BP$Count,ego_result_CC$Count,ego_result_MF$Count),
  Description = c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  type = c(rep("biological process", display_number[1]),
           rep("cellular component", display_number[2]),
           rep("molecular function", display_number[3]))
)

# 5.Select the first 5 words of the path as the name of the path
for (i in 1:nrow(go_enrich_df)){
  description_splite = strsplit(go_enrich_df$Description[i], split = " ")
  description_collapse = paste(description_splite[[1]][1:5], collapse = " ")
  go_enrich_df$Description[i] = description_collapse
  go_enrich_df$Description = gsub(pattern = "NA","",go_enrich_df$Description)
}

# 6.Draw GO bar chart
go_enrich_df$type_order =base::factor(rev(as.integer(rownames(go_enrich_df))), labels = rev(go_enrich_df$Description))

COLS = c('#66C3A5','#8DA1CB','#FD8D62')

ggplot(data = go_enrich_df,aes(x=type_order,y=GeneNumber,fill=type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() +
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched Go Terms") +
  theme_bw()






