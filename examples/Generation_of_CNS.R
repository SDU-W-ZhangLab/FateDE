####################################################################
#########################   SDE model   ###########################
####################################################################
library("diffeqr")

de <- diffeqr::diffeq_setup()

rainbow_cor=rainbow(12)
names(rainbow_cor)=gene_names
#################
# CNS simulation
#################


f <- function(u,p,t) {
  du1=-p[4]*u[1]  #x1---Pax6
  du2=p[2]*u[1]^p[5]/(1+u[1]^p[5]+u[6]^p[9])-p[6]*u[2]           #x2---Mash1
  du3=p[8]*u[2]^p[5]/(1+u[2]^p[5])-p[6]*u[3]                #x3---Zic1
  du4=p[3]*u[3]^p[5]/(1+u[3]^p[5]+u[8]^p[5])-p[6]*u[4]           #x4---Brn2
  du5=p[3]*(u[3]^p[5]+u[4]^p[5])/(1+u[3]^p[5]+u[4]^p[5])-p[6]*u[5]     #x5---Tuj1
  du6=p[2]*u[1]^p[5]/(1+u[1]^p[5]+u[2]^p[9])-p[6]*u[6]             #x6---Hes5
  du7=p[3]*u[6]^p[5]/(1+u[6]^p[5]+u[8]^p[7])-p[6]*u[7]           #x7---Scl
  du8=p[3]*u[6]^p[5]/(1+u[6]^p[5]+u[7]^p[7])-p[6]*u[8]            #x8---Olig2
  du9=p[3]*u[7]^p[5]/(1+u[7]^p[5])-p[6]*u[9]                  #x9---Stat3
  du10=p[3]*u[9]^p[5]/(1+u[9]^p[5])-p[6]*u[10]              #x10---Aldh1L
  du11=p[3]*u[8]^p[5]/(1+u[8]^p[5])-p[6]*u[11]                #x11---Myt1L
  du12=p[3]*u[8]^p[5]/(1+u[8]^p[5])-p[6]*u[12]                 #x12---Sox8

  return(c(du1,du2,du3,du4,du5,du6,du7,du8,du9,du10,du11,du12))
}

g <- function(u,p,t) {
  return(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
}

p <- c(2.2,4,4,0.02,4,1,6,3,8) #p[1]->a_s   p[2]->a  p[3]->a_e  p[4]->η  p[5]->n  p[6]->k
u0 <- c(6,0,0,0,0,0,0,0,0,0,0,0)
tspan <- c(0.0,10.0)


dim1=c("Pax6","Mash1","Zic1","Brn2","Tuj1","Hes5","Scl","Olig2","Stat3","Aldh1L","Myt1L","Sox8")
dim2=paste("cell",1:1001,sep = "_")
dim3=paste("replicate",1:200,sep="_")

simulation=array(NA,dim=c(12,1001,200),dimnames =list(dim1,dim2,dim3) )
for ( i in 1:200){
  print(i)
  prob <- de$SDEProblem(f,g,u0,tspan,p)
  sol <- de$solve(prob,saveat=0.01)
  udf <- as.data.frame(t(sapply(sol$u,identity)))
  simulation[,,i]=t(udf)
}

