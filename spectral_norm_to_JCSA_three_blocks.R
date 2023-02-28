#Simulation for SBM with 3 blocks
#Choosing parameters
n1 <- 1000
p12 <-0.005
p13 <- 0.003
p23 <- 0.002
p1 <- 0.3
p2 <- 0.2
p3 <- 0.1
alpha1 <- 0.3
alpha2 <- 0.5
alpha3 <- 0.2
##Compute zeta
psi1 <- (p1*(1-p1))/(p1*(1-p1))
psi2 <- (p2*(1-p2))/(p1*(1-p1))
psi3 <- (p3*(1-p3))/(p1*(1-p1))
psi12 <- p12*(1-p12)/(p1*(1-p1))
psi13 <- p12*(1-p13)/(p1*(1-p1))
psi23 <- p23*(1-p23)/(p1*(1-p1))
psi0 <- max(psi12,psi13,psi23)
psi1*alpha1
psi2*alpha2
psi3*alpha3
#create SBM
gamma <- (1/sqrt(n1*p1*(1-p1)))
B <- matrix(c(p1,p12,p13,p12,p2,p23,p13,p23,p3),3,3,byrow = TRUE)
Theta1.1 <- matrix(rep(c(1,0,0),n1*alpha1),n1*alpha1,3,byrow = TRUE)
Theta1.2 <- matrix(rep(c(0,1,0),n1*alpha2),n1*alpha2,3,byrow = TRUE)
Theta1.3 <- matrix(rep(c(0,0,1),n1*alpha3),n1*alpha3,3,byrow = TRUE)
Theta1 <- rbind(Theta1.1,Theta1.2,Theta1.3)
##4: generate matrix P=Theta*B*Theta^T
#create matrix P
P1 <- Theta1 %*% B %*% t(Theta1)
require(Rlab)
set.seed(88)
tmp_matrix=matrix(rbern(n1*n1,P1),n1,n1)
tmp_matrix2=tmp_matrix
tmp_matrix2[upper.tri(tmp_matrix2)]=t(tmp_matrix)[upper.tri(tmp_matrix)]
require(gplots)
A <- tmp_matrix2
A_hat <- gamma*A
A_bar <- gamma*P1
A_theota <- A_hat - A_bar
# eigenvalues of A_theota and A_hat
e1 <- eigen(A_theota)$values
e2 <- eigen(A_hat)$values
extrempoint <- 2*sqrt(max(alpha1*psi1,alpha2*psi2,alpha3*psi3)+max((1-alpha1)*psi0,(1-alpha2)*psi0,(1-alpha3)*psi0))
###############################################################################################################################
#histogram of eigenvalues of A_theota
hist(e1,xlim=c(-1.5,1.5),ylim = c(0,35), xlab= NULL, ylab = "Frequency", probability = FALSE,breaks = 40,main=NULL)# 
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)


#histogram of eigenvalues of A_hat
hist(e1,xlim=c(-1.5,8),ylim = c(0,35), xlab= NULL, ylab = "Frequency", probability = FALSE,breaks = 40,main=NULL)# 
#lines(x_var2,y_var2,add=TRUE,col ="blue",lwd=2)
points(e2[1],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[2],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[3],0,col="black", lwd=10,pch=16,add=TRUE)
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
##################################################################################
#Run 1000 simulations

A_hat <- function(p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3){
  n1
  p1
  p2
  p3
  p12
  p13
  p23
  alpha1
  alpha2
  alpha3
  gamma <- (1/sqrt(n1*p1*(1-p1)))
  B <- matrix(c(p1,p12,p13,p12,p2,p23,p13,p23,p3),3,3,byrow = TRUE)
  ##3: Create matrix Theta
  Theta1.1 <- matrix(rep(c(1,0,0),n1*alpha1),n1*alpha1,3,byrow = TRUE)
  Theta1.2 <- matrix(rep(c(0,1,0),n1*alpha2),n1*alpha2,3,byrow = TRUE)
  Theta1.3 <- matrix(rep(c(0,0,1),n1*alpha3),n1*alpha3,3,byrow = TRUE)
  dim(Theta1.1)
  Theta1 <- rbind(Theta1.1,Theta1.2,Theta1.3)
  P1 <- Theta1 %*% B %*% t(Theta1)
  require(Rlab)
  tmp_matrix=matrix(rbern(n1*n1,P1),n1,n1)
  tmp_matrix2=tmp_matrix
  tmp_matrix2[upper.tri(tmp_matrix2)]=t(tmp_matrix)[upper.tri(tmp_matrix)]
  require(gplots)
  A <- tmp_matrix2
  A_hat <- gamma*A
  return(A_hat)
}

#eigenvalue 3
list_eigenvalue3_SBM3 <- c()
for (i in 1:1000) {
  list_eigenvalue3_SBM3[i] <- eigen(A_hat(p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3),symmetric = TRUE)$values[3]
}

list_eigenvalue3_SBM3

sum(list_eigenvalue3_SBM3 >extrempoint)

##########################
#eigenvalue 4
list_eigenvalue4_SBM3 <- c()
for (i in 1:1000) {
  list_eigenvalue4_SBM3[i] <- eigen(A_hat(p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3),symmetric = TRUE)$values[4]
}

list_eigenvalue4_SBM3

#list_eigenvalue_SBM2[5]

sum(list_eigenvalue4_SBM3 >extrempoint)
write.csv(list_eigenvalue4_SBM3, "list_eigenvalue4_SBM3")
write.csv(list_eigenvalue3_SBM3, "list_eigenvalue3_SBM3")
