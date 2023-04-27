library(pheatmap)
library("reticulate")
pd <- import("pandas")
source_python("base.py")


##' Fitting gene expression of one branch
##'
##' @param df1_notfit Branch 1 unfitted data
##'
##' @return a dataframe of fitting
##' @export

fateDE_fit_one = function(df1_notfit){

  backup_expression_family <- VGAM::negbinomial
  fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)"

  # Integerization
  df1_integer = as.data.frame(apply(df1_notfit[,3:ncol(df1_notfit)],1:2,round))
  df1_integer = cbind(df1_notfit$Pseudotime,df1_integer)

  pse_off = seq(min(df1_integer['Pseudotime']),max(df1_integer['Pseudotime']),length.out=nrow(df1_integer))
  pse_off = data.frame(Pseudotime=pse_off)
  df1_fit = data.frame(Pseudotime=pse_off, State = df1_notfit$State)

  print("error genes:")
  for (gene in colnames(df1_notfit)[3:ncol(df1_notfit)]){
    modelFormulaStr <- paste(gene, fullModelFormulaStr,
                             sep = "")
    tryCatch(
      {

        FM_fit=VGAM::vglm(stats::as.formula(modelFormulaStr),
                          family = backup_expression_family,
                          epsilon = 2.0,
                          checkwz = TRUE,
                          data = df1_integer)

        a=monocle::responseMatrix(list(FM_fit),pse_off,response_type="response")
        df1_fit[,gene] = a[1,]


      }, error=function(e){

        print(gene)
      }
    )

  }

  return(df1_fit)

}


##' Identify index of pre somatic cell cells
##'
##' @param exprs_mat Gene expression matrix, cells by now, genes by column
##'
##' @return Index of pre somatic cell cells
##' @export

Determine_transition_state=function(exprs_mat){
  ICs=c()
  cut_index=cut(1:nrow(exprs_mat),breaks =100 ,include.lowest = T)
  index=unique(cut_index)
  for(i in 1:length(index)){
    exprs_mat_i=exprs_mat[which(cut_index==index[i]),]
    cor_genes=stats::cor(exprs_mat_i)
    cor_cells=stats::cor(t(exprs_mat_i))
    ic=mean(cor_genes[lower.tri(cor_genes)])/mean(cor_cells[lower.tri(cor_cells)])
    ICs=c(ICs,ic)
  }
  plot(ICs)
  transition_index=nrow(exprs_mat[sapply(cut_index,function(x) x%in% index[1:which.max(ICs)]),])

  return(transition_index)
}



##' FateDE: Searching for Destiny Determinants
##'
##' @param branch1 A Dataframe,cells by row, the first column is Pseudotime, the second column is cell State,
##' other columns are genes.
##' @param branch2 The same as branch1.
##' @param G1 A vector,the set of genes that increase in expression over pseudotime in branch 1.
##' @param G2 A vector,the set of genes that increase in expression over pseudotime in branch 2.
##' @param GP_FSV A number, screening threshold in Gaussian process.
##' @param logN A logic value, is a logified expression matrix used in Gaussian processes.
##' @param TFs A vector, extract the TF of the expression matrix for calculation.
##' @param cluster_meth clustering method, clustering methods using the same as pheatmap.
##'
##' @return A list, include G1, G2, GP's results, G1_1, G2_1, one_two_cor_pairwise_all's results
##' @export

FateDE_main=function(branch1, branch2, G1=NULL, G2=NULL,
                     GP_FSV = 0.9,
                     logN=T,
                     TFs = NULL,
                     cluster_meth = "ward.D2"){


  if(is.null(branch1) || is.null(branch2)){
    stop("Please enter the branch presentation matrix.")
  }

  branch1 = branch1[order(branch1$Pseudotime),]
  branch2 = branch2[order(branch2$Pseudotime),]
  ################################ step1: Obtain a constantly increasing set of genes ################################
  #############################################################################################
  if(is.null(G1) && is.null(G2)){

    cat("############# Obtain genes with consistently increasing expression in branch 1 #############\n")
    cat("\n############################################################################################\n")
    # Cluster Branch 1
    #################
    htmap<- pheatmap::pheatmap(t(branch1[,3:ncol(branch1)]),
                     cluster_cols = FALSE,
                     clustering_method = cluster_meth,
                     treeheight_row=100,
                     show_colnames=F)
    # Cluster Tree for Extracting Row Directions (Genes) from Heat Maps
    cluster <- htmap$tree_row
    # Cluster the clustering tree
    cluster_num<-readline(prompt = "How many clusters do you want to divide: ")
    cut <- stats::cutree(cluster,cluster_num)
    annotation_row = data.frame(cluster =as.character(cut))
    rownames(annotation_row) = names(cut)
    htmap<- pheatmap::pheatmap(t(branch1[,3:ncol(branch1)]),
                     cluster_cols = FALSE,
                     clustering_method = cluster_meth,
                     treeheight_row=100,
                     show_colnames=F,
                     annotation_row = annotation_row)
    # Obtain branch 1 genes that increase
    G1 = c()
    cluster_if = "T"
    while(cluster_if == "T"){

      cluster_if<-readline(prompt = "Do you need to continue taking cluster,T/F: ")
      if(cluster_if == "F"){
        break
      }
      the_cluster_num<-readline(prompt = "Which cluster do you want to obtain: ")
      G1 = c(G1, names(cut)[cut==the_cluster_num])
    }


    cat("\n############################################################################################\n")
    cat("\n############# Obtain genes with consistently increasing expression in branch 2 #############\n")
    cat("\n############################################################################################\n")

    # Cluster Branch 1
    #################
    htmap<- pheatmap::pheatmap(t(branch2[,3:ncol(branch2)]),
                     cluster_cols = FALSE,
                     clustering_method = cluster_meth,
                     treeheight_row=100,
                     show_colnames=F)
    # Cluster Tree for Extracting Row Directions (Genes) from Heat Maps
    cluster <- htmap$tree_row
    # Cluster the clustering tree
    cluster_num<-readline(prompt = "How many clusters do you want to divide: ")
    cut <- stats::cutree(cluster,cluster_num)
    annotation_row = data.frame(cluster =as.character(cut))
    rownames(annotation_row) = names(cut)
    htmap<- pheatmap::pheatmap(t(branch2[,3:ncol(branch2)]),
                     cluster_cols = FALSE,
                     clustering_method = cluster_meth,
                     treeheight_row=100,
                     show_colnames=F,
                     annotation_row = annotation_row)
    # Obtain branch 1 genes that increase
    G2 = c()
    cluster_if = "T"
    while(cluster_if == "T"){

      cluster_if<-readline(prompt = "Do you need to continue taking cluster,T/F: ")
      if(cluster_if == "F"){
        break
      }
      the_cluster_num<-readline(prompt = "Which cluster do you want to obtain: ")
      G2 = c(G2, names(cut)[cut==the_cluster_num])
    }

  } else{
    G1 = G1
    G2 = G2
  }

  ################ step2: Obtain a set of genes that increase first and then decrease later ################################
  #############################################################################################

  G = find_genes(branch1, branch2, intersect(G1,colnames(branch2)),
                 intersect(G2,colnames(branch1)), "fate1", "fate2", logN=logN)
  names(G) = c('G1_GP','G2_GP')
  G1_1 = G[[1]]
  G2_1 = G[[2]]

  G1_1 = subset(G1_1, G1_1$l < 3.0)
  G1_1 = subset(G1_1, G1_1$FSV > GP_FSV)
  G1_1 = c(G1_1$gene)

  if(is.null(G1_1)){
    stop("G1_1 is null, please lower GP_FSV.")
  }

  for (gene in G1_1) {  # Remove the constant increase in df2
    a = branch2[[gene]]
    gene_max_index = which(a==max(a), arr.ind = TRUE)
    gene_max_pse = branch2[gene_max_index,'Pseudotime']
    if(gene_max_pse == max(branch2[['Pseudotime']])){
      # print('Continuously decreasing')
      # print(gene)
      G1_1 = G1_1[-which(G1_1 == gene)]
    }
  }


  G2_1 = subset(G2_1, G2_1$l < 3.0)
  G2_1 = subset(G2_1, G2_1$FSV > GP_FSV)
  G2_1 = c(G2_1$gene)

  if(is.null(G2_1)){
    stop("G2_1 is null, please lower GP_FSV.")
  }

  for (gene in G2_1) {  # Remove the constant increase in df1
    a = branch1[[gene]]
    gene_max_index = which(a==max(a), arr.ind = TRUE)
    gene_max_pse = branch1[gene_max_index,'Pseudotime']
    if(gene_max_pse == max(branch1[['Pseudotime']])){
      # print('Continuously decreasing')
      # print(gene)
      G2_1 = G2_1[-which(G2_1 == gene)]
    }
  }


  ################################ step3: Correlation coefficient analysis ################################
  #############################################################################################
  if(!is.null(TFs)){
    G1_1 = intersect(G1_1,TFs) # Increase in df1 and then decrease in df2
    G2_1 = intersect(G2_1,TFs) # Increase in df2 and then decrease in df1
  }

  gene_pair_set=merge(G1_1,G2_1) # Take Cartesian product

  progenitors1 = subset(branch1, branch1$State==unique(branch1$State)[1]) # Ensure cell state
  differentiations1 = subset(branch1, branch1$State==unique(branch1$State)[2])

  progenitors2 = subset(branch2, branch2$State==unique(branch2$State)[1])
  differentiations2 = subset(branch2, branch2$State==unique(branch2$State)[2])

  one_two_cor_pairwise_all=c()
  for( i in 1:nrow(gene_pair_set)){

    gene_pair=unlist(gene_pair_set[i,])
    one_two_cor_pairwise=get_one_two_cor_pairwise(branch1,branch2,progenitors1,differentiations1,progenitors2,differentiations2,gene_pair)
    one_two_cor_pairwise_all=data.frame(rbind(one_two_cor_pairwise_all,one_two_cor_pairwise))
  }

  rownames(one_two_cor_pairwise_all)=apply(gene_pair_set,1, function (x) paste(x,collapse  = "_"))
  colnames(one_two_cor_pairwise_all)=c("branch1_observation","branch2_observation","one_order_progenitor","one_order_differentiations","two_order_branch1","two_order_branch2")

  result = list(G1, G2, G, G1_1, G2_1, as.data.frame(one_two_cor_pairwise_all))
  names(result) = c("G1", "G2", "GP", "G1_1", "G2_1", "one_two_cor_pairwise_all")
  return(result)
}

##' Calculating Score
##'
##' @param branch1 branch1 cells by row, the first column is Pseudotime, the second column is cell State,
##' other columns are genes.
##' @param branch2 the same as branch1
##' @param progenitors1 Precursor somatic cell cell of branch 1
##' @param differentiations1 Differentiated cells of branch 1
##' @param progenitors2 Precursor somatic cell cell of branch 2
##' @param differentiations2 Differentiated cells of branch 2
##' @param gene_pair A pair of genes to be calculated
##'
##' @return Score for gene pair

get_one_two_cor_pairwise=function(branch1,branch2,progenitors1,differentiations1,progenitors2,differentiations2,gene_pair){

  branch1_genepair=branch1[,gene_pair]
  branch2_genepair=branch2[,gene_pair]

  progenitors1_genepair=progenitors1[,gene_pair]
  progenitors2_genepair=progenitors2[,gene_pair]

  differentiations1_genepair=differentiations1[,gene_pair]
  differentiations2_genepair=differentiations2[,gene_pair]

  #control: correlation of exp in two genes at each branch
  one_order_branch1=abs(stats::cor(branch1_genepair)[1,2])
  one_order_branch2=abs(stats::cor(branch2_genepair)[1,2])


  one_order_cor_progenitors=stats::cor(rbind(progenitors1_genepair,progenitors2_genepair))[1,2]
  one_order_cor_differentiations=stats::cor(rbind(differentiations1_genepair,differentiations2_genepair))[1,2]

  branch1_transition=nrow(progenitors1)
  branch2_transition=nrow(progenitors2)

  ##  time series of diff in branch1
  diff_1=branch1_genepair[-1,]-branch1_genepair[-nrow(branch1_genepair),]
  diff_2=branch2_genepair[-1,]-branch2_genepair[-nrow(branch2_genepair),]

  two_order_cor_branch1=stats::cor(Fit_t(abs(diff_1[,gene_pair[1]])),Fit_t(abs(diff_1[,gene_pair[2]])))
  two_order_cor_branch2=stats::cor(Fit_t(abs(diff_2[,gene_pair[1]])),Fit_t(abs(diff_2[,gene_pair[2]])))

  one_two_cor_pairwise=c(one_order_branch1,one_order_branch2,one_order_cor_progenitors,one_order_cor_differentiations,two_order_cor_branch1,two_order_cor_branch2)
  names(one_two_cor_pairwise)=c("one_order_branch1","one_order_branch2","one_order_progenitor","one_order_differentiations","two_order_branch1","two_order_branch2")

  return(one_two_cor_pairwise)


}


Fit_t=function(exp_value){

  d=data.frame(cbind(1:length(exp_value),exp_value))
  colnames(d)=c("x","y")
  #tmp <- lm(formula = y ~ x+1 )
  #tmp<- nls(y ~ k*(x^z)*exp(-x),data.frame(d),weights = c(rep(100,10),rep(1,length(exp_value)-10)),start = list(k = 4, z=2 ))
  #coeffi=coef(tmp)
  #exp_fitted=predict(tmp,list(x=d[,1]))
  tmp <- stats::lm(y ~ poly(x,6),data=d)
  exp_fitted=tmp$fitted.value
  return(exp_fitted)

}



