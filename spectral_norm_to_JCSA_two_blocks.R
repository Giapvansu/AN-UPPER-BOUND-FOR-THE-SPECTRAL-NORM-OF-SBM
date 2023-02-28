#SBM with two blocks
#Choosing parameters
n1 <- 1000
p0 <- 0.02
p1 <- 0.3
p2 <- 0.1
rate <- 0.45
##Compute zeta0 and zeta2
zeta1 <- p1*(1-p1)/(p1*(1-p1))
zeta0 <- (p0*(1-p0))/(p1*(1-p1))
zeta2 <- (p2*(1-p2))/(p1*(1-p1))
#Create the Stochastic block matrix P 
n2 <- rate*n1
n3 <- n1-n2
gamma <- 1/sqrt(n1*p1*(1-p1))
B <- matrix(c(p1,p0,p0,p2),2,2,byrow = TRUE)
##3: Create matrix Theta
Theta1.1 <- matrix(rep(c(1,0),n2),n2,2,byrow = TRUE)
Theta1.2 <- matrix(rep(c(0,1),n3),n3,2,byrow = TRUE)
dim(Theta1.2)
Theta1 <- rbind(Theta1.1,Theta1.2)
##4: generate matrix P=Theta*B*Theta^T
#create matrix P
P1 <- Theta1 %*% B %*% t(Theta1)
require(Rlab)
#set.seed(88)
tmp_matrix=matrix(rbern(n1*n1,P1),n1,n1)
tmp_matrix2=tmp_matrix
tmp_matrix2[upper.tri(tmp_matrix2)]=t(tmp_matrix)[upper.tri(tmp_matrix)]
require(gplots)
A <- tmp_matrix2
A_hat <- gamma*A
A_bar <- gamma*P1
A_theota <- A_hat - A_bar
###############################################################################################
#compute upper bound
extrempoint <- 2*sqrt(max(rate*zeta1,rate*zeta2)+max((1-rate)*zeta0,rate*zeta0))
#Compute eigenvalues set of matrix A_theota and A_hat
e <- eigen(A_theota)$values
e2 <- eigen(A_hat)$values
#Histogram of eigenvalues of A_theota
hist(e,xlim=c(-1.5,1.5),ylim=c(0,60),xlab= NULL,,ylab = "Frequency",probability = FALSE,breaks = 30,main = NULL) #main = paste("p1=",p1,",","p2=", p2,",p0= 2log(N)/N^2,","alpha=", rate ))
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
#Histogram of eigenvalues of A_hat
hist(e,xlim=c(-1.5,10),ylim = c(0,60), xlab= NULL, ylab = "Frequency", probability = FALSE,breaks = 30,main=NULL)# main = paste("p1=", p1_1,"," ,"p2=",p2_2 ,"p3=",p3_2 ,",","p0=log(N)/N"))
points(e2[1],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[2],0,col="black", lwd=10,pch=16,add=TRUE)
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
heatmap(A)
heatmap(A_theota)
########################################################
#Run 1000 simulations

A_hat <- function(n1,p1,p2,p0,rate){
  n1
  p1
  p2
  p0
  rate
  n2 <- rate*n1
  n3 <- n1-n2
  B <- matrix(c(p1,p0,p0,p2),2,2,byrow = TRUE)
  Theta1.1 <- matrix(rep(c(1,0),n2),n2,2,byrow = TRUE)
  Theta1.2 <- matrix(rep(c(0,1),n3),n3,2,byrow = TRUE)
  dim(Theta1.2)
  Theta1 <- rbind(Theta1.1,Theta1.2)
  P1 <- Theta1 %*% B %*% t(Theta1)
  require(Rlab)
  tmp_matrix=matrix(rbern(n1*n1,P1),n1,n1)
  tmp_matrix2=tmp_matrix
  tmp_matrix2[upper.tri(tmp_matrix2)]=t(tmp_matrix)[upper.tri(tmp_matrix)]
  A <- tmp_matrix2
  gamma <- 1/sqrt(n1*p1*(1-p1))
  A_hat <- gamma*A
  return(A_hat)
}

#list of second eigenvalues of A_hat
list_eigenvalue2_SBM2 <- c()
for (i in 1:1000) {
  list_eigenvalue2_SBM2[i] <- eigen(A_hat(n1,p1,p2,p0,rate),symmetric = TRUE)$values[2]
}

list_eigenvalue2_SBM2
# number of second eigenvalues outside upper bound
sum(list_eigenvalue2_SBM2 >extrempoint)
#list of third eigenvalues of A_hat
#eigenvalue 3
list_eigenvalue3_SBM2 <- c()
for (i in 1:1000) {
  list_eigenvalue3_SBM2[i] <- eigen(A_hat(n1,p1,p2,p0,rate),symmetric = TRUE)$values[3]
}

list_eigenvalue1_SBM2

sum(list_eigenvalue3_SBM2 >extrempoint)
write.csv(list_eigenvalue2_SBM2, "list_eigenvalue2_SBM2")
write.csv(list_eigenvalue3_SBM2, "list_eigenvalue3_SBM2")
