library(FateDE)

data(simulation)

#obtain the corresponding lineage for each simulation run:
neuron_cell_ids <- which(simulation['Mash1', 1001, ] > 3)
astrocyte_cell_ids <- which(simulation['Scl',1001 , ] > 3)
oligodendrocyte_cell_ids <- which(simulation['Olig2', 1001, ] > 3)

gene_names <- c("Pax6","Mash1","Zic1","Brn2","Tuj1","Hes5","Scl","Olig2","Stat3","Aldh1L","Myt1L","Sox8")
neuron_exprs_mat <- as.data.frame(t((simulation[, , sample(neuron_cell_ids,1)])))
astrocyte_exprs_mat <-  as.data.frame(t(simulation[, , sample(astrocyte_cell_ids,1)]))
oligodendrocyte_exprs_mat <-  as.data.frame(t(simulation[, , sample(oligodendrocyte_cell_ids,1)]))


neuron_trans=Determine_transition_state(neuron_exprs_mat)
Pseudotime = seq(1,nrow(neuron_exprs_mat),by=1)
State = c(rep("State 1",neuron_trans), rep("State 2",nrow(neuron_exprs_mat)-neuron_trans))
neuron_exprs_mat = cbind(Pseudotime,State,neuron_exprs_mat)

astrocyte_trans=Determine_transition_state(astrocyte_exprs_mat)
Pseudotime = seq(1,nrow(astrocyte_exprs_mat),by=1)
State = c(rep("State 1",astrocyte_trans), rep("State 2",nrow(astrocyte_exprs_mat)-astrocyte_trans))
astrocyte_exprs_mat = cbind(Pseudotime,State,astrocyte_exprs_mat)

oligodendrocyte_trans=Determine_transition_state(oligodendrocyte_exprs_mat)
Pseudotime = seq(1,nrow(oligodendrocyte_exprs_mat),by=1)
State = c(rep("State 1",oligodendrocyte_trans), rep("State 2",nrow(oligodendrocyte_exprs_mat)-oligodendrocyte_trans))
oligodendrocyte_exprs_mat = cbind(Pseudotime,State,oligodendrocyte_exprs_mat)

# N-A
result_NA = FateDE_main(neuron_exprs_mat, astrocyte_exprs_mat, logN=F, GP_FSV = 0.3)
# N-O
result_NO = FateDE_main(neuron_exprs_mat, oligodendrocyte_exprs_mat, logN=F, GP_FSV = 0.3)
# A-O
result_AO = FateDE_main(astrocyte_exprs_mat, oligodendrocyte_exprs_mat, logN=F, GP_FSV = 0.3)

##############################
#figure1: time series raw data
##############################
par(mfrow=c(3,2),mar=c(4,4,4,4))
gene_pair = c("Mash1","Hes5")
##raw time series and fitted time series in branch1
plot(1,type="n",xlim=c(0,1000),ylim=c(-0.05,5),xlab="",ylab="",frame.plot=F,main="Branch1",las=1)
points(neuron_exprs_mat[,gene_pair[1]],type="l",lwd=2,lty=1, col="darkred")
#plot(1,type="n",xlim=c(0,1000),ylim=c(-0.05,2.5),xlab="",ylab="",frame.plot=F,main="gene2")
points(neuron_exprs_mat[,gene_pair[2]],type="l",lwd=2,lty=4, col="darkred")
legend(550,3,c(gene_pair[1],gene_pair[2]),col=c("darkred","darkred"),lty=c(1,4),cex=0.6,bty="n")
abline(v=neuron_trans,lty=2)

##raw time series and fitted time series in branch2
plot(1,type="n",xlim=c(0,1000),ylim=c(-0.05,5),xlab="",ylab="",frame.plot=F,main="Branch2",las=1)
points(astrocyte_exprs_mat[,gene_pair[1]],type="l",lwd=2,lty=1, col="darkblue")
points(astrocyte_exprs_mat[,gene_pair[2]],type="l",lwd=2,lty=4, col="darkblue")
legend(550,3,c(gene_pair[1],gene_pair[2]),col=c("darkblue","darkblue"),lty=c(1,4),cex=0.6,bty="n")
abline(v=astrocyte_trans,lty=2)

#par(mfrow=c(1,1))
###############################
# figure3:  pairwise scatter plot
###############################

##raw gene pair in progenitors

plot(1,type="n",xlim=c(0,5),ylim=c(0,5),frame.plot=F,xlab=gene_pair[1],ylab=gene_pair[2],las=1) ##

points(neuron_exprs_mat[1:neuron_trans ,gene_pair],pch=16,cex=1.5,col="grey")
points(astrocyte_exprs_mat[1:astrocyte_trans ,gene_pair],pch=16,cex=1.5,col="grey")

#plot(1,type="n",xlim=c(0,2.5),ylim=c(0,2.5),frame.plot=F,xlab=gene_pair[1],ylab=gene_pair[2])

points(neuron_exprs_mat[neuron_trans:nrow(neuron_exprs_mat) ,gene_pair],pch=16,cex=1.5,col="darkred")
points(astrocyte_exprs_mat[astrocyte_trans:nrow(neuron_exprs_mat) ,gene_pair],pch=16,cex=1.5,col="darkblue")
text(0.6,3,paste("r1=",round(result_NA$one_two_cor_pairwise_all["Mash1_Hes5","one_order_progenitor"],3)),cex=0.5)
text(0.6,2.5,paste("r2=",round(result_NA$one_two_cor_pairwise_all["Mash1_Hes5","one_order_differentiations"],3)),cex=0.5) #####
legend(3,3.5,c("progenitors","branch1","branch2"),col=c("grey","darkred","darkblue"),lty=1,bty="n",cex=0.6)

###############################
# figure2 plot pairwise diff
###############################
branch1_genepair = neuron_exprs_mat[,gene_pair]
branch2_genepair = astrocyte_exprs_mat[,gene_pair]

##  time series of diff in branch1
diff_1=branch1_genepair[-1,]-branch1_genepair[-nrow(branch1_genepair),]
diff_2=branch2_genepair[-1,]-branch2_genepair[-nrow(branch2_genepair),]

two_order_cor_branch1=cor(Fit_t(abs(diff_1[,gene_pair[1]])),Fit_t(abs(diff_1[,gene_pair[2]])))
two_order_cor_branch2=cor(Fit_t(abs(diff_2[,gene_pair[1]])),Fit_t(abs(diff_2[,gene_pair[2]])))


plot(1,type="n",xlim=c(0,1000),ylim=c(0,0.04),xlab="",ylab="",frame.plot=F,las=1)
#points(abs(diff_1[,gene_pair[1]]),type="l",lwd=1,lty=1, col= brewer.pal(8, 'Reds')[3])
points(Fit_t(abs(diff_1[,gene_pair[1]])),type="l",lwd=4,lty=1,col= "darkred")


#plot(1,type="n",xlim=c(0,1000),ylim=c(0,0.05),xlab="",ylab="",frame.plot=F,main=gene_pair[2],las=1)
#points(abs(diff_1[,gene_pair[2]]),type="l",lwd=1,lty=1, col= brewer.pal(8, 'Reds')[3])
points(Fit_t(abs(diff_1[,gene_pair[2]])),type="l",lwd=4,lty=4,col= "darkred")
text(500,0.035,paste("r3=",round(two_order_cor_branch1,3)),cex=0.8)
legend(700,0.035,c(gene_pair[1],gene_pair[2]),col=c("darkred","darkred"),lty=c(1,4),cex=0.6,bty="n")


##time series of diff in branch2


plot(1,type="n",xlim=c(0,1000),ylim=c(0,0.04),xlab="",ylab="",frame.plot=F,las=1)
#points(abs(diff_2[,gene_pair[1]]),type="l",lwd=1,lty=1, col= brewer.pal(8, 'Blues')[4])
points(Fit_t(abs(diff_2[,gene_pair[1]])),type="l",lwd=4,lty=1,col= "darkblue")


#plot(1,type="n",xlim=c(0,1000),ylim=c(0,0.05),xlab="",ylab="",frame.plot=F,main=gene_pair[2],las=1)
#points(abs(diff_2[,gene_pair[2]]),type="l",lwd=1,lty=1, col= brewer.pal(8, 'Blues')[4])
points(Fit_t(abs(diff_2[,gene_pair[2]])),type="l",lwd=4,lty=4,col= "darkblue")
text(500,0.035,paste("r4=",round(two_order_cor_branch2,3)),cex=0.8)
legend(700,0.035,c(gene_pair[1],gene_pair[2]),col=c("darkblue","darkblue"),lty=c(1,4),cex=0.6,bty="n")



Fit_t=function(exp_value){

  d=data.frame(cbind(1:length(exp_value),exp_value))
  colnames(d)=c("x","y")
  #tmp <- lm(formula = y ~ x+1 )
  #tmp<- nls(y ~ k*(x^z)*exp(-x),data.frame(d),weights = c(rep(100,10),rep(1,length(exp_value)-10)),start = list(k = 4, z=2 ))
  #coeffi=coef(tmp)
  #exp_fitted=predict(tmp,list(x=d[,1]))
  tmp <- lm(y ~ poly(x,6),data=d)
  exp_fitted=tmp$fitted.value
  return(exp_fitted)

}


###############################
# figure plot circular diagram
###############################
library('ComplexHeatmap')

library('circlize')

library("RColorBrewer")

library(dendextend)

#转化matrix格式矩阵及数据归一化
cir1 <- as.matrix(abs(astrocyte_oligodendrocyt[,3:6]))

cir2 <- as.matrix(abs(neuron_astrocyt[,3:6]))

cir3 <- as.matrix(abs(neuron_oligodendrocyt[,3:6]))


#颜色设定：
mycol=colorRamp2(c(0.5,0.9,1),c("blue","white","red"))
circos.clear()

#cir1分组热图绘制#

ann_row = data.frame(pathway=c(rep("pathway1",3),rep("pathway2",3),rep("pathway3",3)))#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir1)

ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix

#分组绘图

circos.par(gap.after=c(2,2,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir1,col=mycol,

               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈

               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外

               track.height = 0.2, #轨道的高度，数值越大圆环越粗

               rownames.col="black",

               bg.border="black", #背景边缘颜色

               split = ann_row,#用行注释分裂热图

               show.sector.labels = F,

               rownames.cex=0.8,#字体大小

               rownames.font=0.5,#字体粗细

               cluster=FALSE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

               dend.track.height=0.18,#调整行聚类树的高度

               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称

                 color_branches(dend,k=3,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色

               }

)
lg=Legend(title="Legend",col_fun=mycol,direction = c("horizontal"))  #horizontal vertical
grid.draw(lg)

circos.clear()
mycol=colorRamp2(c(0.5,0.9,1),c("blue","white","red"))
#cir2分组热图绘制#

ann_row = data.frame(pathway=c(rep("pathway1",4),rep("pathway2",4),rep("pathway3",4),rep("pathway4",4)))#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir2)

ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix

#分组绘图

circos.par(gap.after=c(2,2,2,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir2,col=mycol,

               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈

               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外

               track.height = 0.2, #轨道的高度，数值越大圆环越粗

               rownames.col="black",

               bg.border="black", #背景边缘颜色

               split = ann_row,#用行注释分裂热图

               show.sector.labels = F,

               rownames.cex=0.8,#字体大小

               rownames.font=0.5,#字体粗细

               cluster=FALSE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

               dend.track.height=0.18,#调整行聚类树的高度

               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称

                 color_branches(dend,k=4,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色

               }

)

lg=Legend(title="Legend",col_fun=mycol,direction = c("horizontal"))  #horizontal vertical
grid.draw(lg)






mycol=colorRamp2(c(0.5,0.9,1),c("blue","white","red"))
circos.clear()

#cir3分组热图绘制#

ann_row = data.frame(pathway=c(rep("pathway1",4),rep("pathway2",4),rep("pathway3",4),rep("pathway4",4)))#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir3)

ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix

#分组绘图

circos.par(gap.after=c(2,2,2,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir3,col=mycol,

               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈

               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外

               track.height = 0.2, #轨道的高度，数值越大圆环越粗

               rownames.col="black",

               bg.border="black", #背景边缘颜色

               split = ann_row,#用行注释分裂热图

               show.sector.labels = F,

               rownames.cex=0.8,#字体大小

               rownames.font=0.5,#字体粗细

               cluster=FALSE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

               dend.track.height=0.18,#调整行聚类树的高度

               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称

                 color_branches(dend,k=3,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色

               }

)

lg=Legend(title="Legend",col_fun=mycol,direction = c("horizontal"))  #horizontal vertical
grid.draw(lg)





