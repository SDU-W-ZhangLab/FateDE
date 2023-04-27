library(igraph)
library(mgcv)
library(quadprog)
library(Rcpp)
library(inline)
library(RcppArmadillo)
library(pbapply)
library(glmnet)
library(parallel)
library(dplyr)
library(plot3D)
library(monocle)
library(crestree)
library(FateDE)

data(crest)
set.seed(123)
nc.cells <- crest$nc.cells
emb <- crest$emb
fpm <- crest$fpm
wgm <- crest$wgm
wgwm <- crest$wgwm


#### Trajectory inference ####
metrics <- "cosine"
M <- length(nc.cells)
lambda <- 150
sigma <- 0.015
z <- ppt.tree(X=fpm[rownames(wgm),nc.cells], emb=emb, lambda=lambda, sigma=sigma, metrics=metrics, M=M, err.cut = 5e-3, n.steps=50, seed=1, plot=FALSE)
plotppt(z,emb,tips=FALSE,cex.tree = 0.1,cex.main=0.2,lwd.tree = 1)

plotppt(z,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)
ppt <- cleanup.branches(z,tips.remove = c(139,295))
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)
ppt <- setroot(ppt,root=205) # 选择root

visualise.trajectory()
#### Cellular mapping ####
cell <- nc.cells[2] # choose a cell
cell = 'p87_O9'
pprobs <- ppt$R[cell,] # probabilities of tree projections
plotppt(ppt,emb,pattern.tree = ppt$R[cell,],cex.tree = 1,lwd.tree = 0.1) # plot probabilities using pattern.tree parameter
points(emb[cell,1],emb[cell,2],cex=1,pch=19,col="black") # show cell position on embedding
ppt <- project.cells.onto.ppt(ppt,emb,n.mapping = 100)

#### Tree dependent genes ####
ppt <- test.associated.genes(ppt,n.map=1,fpm,summary=TRUE,n.cores=1)
head(ppt$stat.association[order(ppt$stat.association$pval),])
genes.tree <- crest$genes.tree

ppt$stat.association$sign <- FALSE
ppt$stat.association[genes.tree,]$sign <- TRUE

ppt <- fit.associated.genes(ppt,fpm,n.map=1,n.cores=1)

gene = 'Neurog2'
visualise.trajectory(ppt,gene,fpm[gene,],cex.main = 3,lwd.t2=0.5)

#### Up-genes in branch 1 ####
gene = 'Neurog2'
# 分支一
plotppt(ppt,emb[,],tips=TRUE,tree.col = ppt$pp.info$color,forks=TRUE,cex.tree = 1,lwd.tree = 0.1) # visualize tree tips
zseg <- extract.subtree(ppt,c("205","304")) # select root and terminal leave of the trajectory
plotppt(ppt,emb,gene=gene,mat=fpm,cex.main=1,cex.tree = 1.5,lwd.tree = 0.1,subtree=zseg)

stat.subtree <- test.associated.genes(ppt,n.map=1,fpm,subtree = zseg,n.cores=1)
cells.subtree <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% zseg$segs]

genes.subtree <- rownames(stat.subtree)[stat.subtree$sign==TRUE & stat.subtree$A > 2]

hc <- hclust(dist(ppt$fit.summary[intersect(rownames(ppt$fit.summary),genes.subtree),cells.subtree]),method="ward.D") # hierarchical clustering
clust <- cutree(hc,4) # partition of genes in 4 clusters

visualise.clusters(ppt,emb,clust=clust,cex.gene=1,cex.cell=0.05,cex.tree=0.2,subtree=zseg)

visualise.trajectory(ppt,"Neurog2",fpm["Neurog2",],cex.main = 3,lwd.t2=0.5,subtree=zseg)

G1=names(clust[clust %in% c(2,4)])  # fate1_logis_genes_ls

#### Up-genes in branch 2 ####
plotppt(ppt,emb[,],tips=TRUE,tree.col = ppt$pp.info$color,forks=TRUE,cex.tree = 1,lwd.tree = 0.1) # visualize tree tips
zseg <- extract.subtree(ppt,c("205","166")) # select root and terminal leave of the trajectory
plotppt(ppt,emb,gene=gene,mat=fpm,cex.main=1,cex.tree = 1.5,lwd.tree = 0.1,subtree=zseg)

stat.subtree <- test.associated.genes(ppt,n.map=1,fpm,subtree = zseg,n.cores=1)
cells.subtree <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% zseg$segs]

genes.subtree <- rownames(stat.subtree)[stat.subtree$sign==TRUE & stat.subtree$A > 2]

hc <- hclust(dist(ppt$fit.summary[intersect(rownames(ppt$fit.summary),genes.subtree),cells.subtree]),method="ward.D") # hierarchical clustering
clust <- cutree(hc,4) # partition of genes in 4 clusters

visualise.clusters(ppt,emb,clust=clust,cex.gene=1,cex.cell=0.05,cex.tree=0.2,subtree=zseg)

visualise.trajectory(ppt,"Twist1",fpm["Twist1",],cex.main = 3,lwd.t2=0.5,subtree=zseg)
G2=names(clust[clust %in% c(2,4)])  # fate2_logis_genes_ls

##### fitting #######
expr = data.frame(t(fpm[, rownames(ppt$cell.summary)]))

tmp = data.frame(
  Pseudotime = ppt$cell.summary$t,
  State = ppt$cell.summary$seg
)

rownames(tmp) = rownames(expr)
tmp = cbind(tmp,expr)

df1_notfit = subset(tmp, milestone_id %in% c(8,10,9,6)) # 分析205-304
df1_notfit = df1_notfit[order(df1_notfit$Pseudotime),]
for (i in 1:nrow(df1_notfit)){

  if (df1_notfit[i,'State'] %in% c(8,10,6)){

    df1_notfit[i,'State'] = 'State 1'

  }  else{
    df1_notfit[i,'State'] = 'State 2'
  }
}

df2_notfit = subset(tmp, milestone_id %in% c(8,2,10,7,5,3,6,1,4)) # 分析205-166
df2_notfit = df2_notfit[order(df2_notfit$Pseudotime),]
for (i in 1:nrow(df2_notfit)){

  if (df2_notfit[i,'State'] %in% c(8,10,6)){

    df2_notfit[i,'State'] = 'State 1'

  }  else{
    df2_notfit[i,'State'] = 'State 3'
  }
}

df1_fit = fateDE_fit_one(df1_notfit)
df2_fit = fateDE_fit_one(df2_notfit)

##### FateDE #######
tfs = read.table('TF/mouse_TF_list.txt',sep = '\n')
result = FateDE_main(df1_fit, df2_fit, G1, G2, GP_FSV = 0.6, TFs =tfs$V1)
G1_1 = result$G1_1
G2_1 = result$G2_1

#################### plot ###########################

# single gene curve
gene = 'Prdm12'
plot_SingleGene_trend(list(df1_notfit,df2_notfit),df1_fit,df2_fit,gene)+
  scale_color_manual(values=c("#999999", "#0073C2FF",'#c74545'))

# double genes curve
double_gene = c('Prdm12','Mef2c')
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


