#Simulation for SBM with 3 blocks
#Choosing parameters
n1 <- 1000
p12 <-0.005
p13 <- 0.005
p14 <- 0.005
p15 <- 0.004
p23 <- 0.004
p24 <- 0.005
p25 <- 0.005
p34 <- 0.003
p35 <- 0.003
p45 <- 0.005
p1 <- 0.25
p2 <- 0.25
p3 <- 0.2
p4 <- 0.15
p5 <- 0.1
alpha1 <- 0.2
alpha2 <- 0.2
alpha3 <- 0.2
alpha4 <- 0.2
alpha5 <- 0.2
##Compute zeta
psi1 <- (p1*(1-p1))/(p1*(1-p1))
psi2 <- (p2*(1-p2))/(p1*(1-p1))
psi3 <- (p3*(1-p3))/(p1*(1-p1))
psi4 <- (p4*(1-p4))/(p1*(1-p1))
psi5 <- (p5*(1-p5))/(p1*(1-p1))
psi12 <- p12*(1-p12)/(p1*(1-p1))
psi13 <- p12*(1-p13)/(p1*(1-p1))
psi14 <- p14*(1-p14)/(p1*(1-p1))
psi15 <- p15*(1-p15)/(p1*(1-p1))
psi23 <- p23*(1-p23)/(p1*(1-p1))
psi24 <- p24*(1-p24)/(p1*(1-p1))
psi25 <- p25*(1-p25)/(p1*(1-p1))
psi34 <- p34*(1-p34)/(p1*(1-p1))
psi35 <- p35*(1-p35)/(p1*(1-p1))
psi45 <- p45*(1-p45)/(p1*(1-p1))
psi0 <- max(psi12,psi13,psi14, psi15,psi23,psi24, psi25,psi34,psi35,psi45)
psi1*alpha1
psi2*alpha2
psi3*alpha3
#create SBM
gamma <- (1/sqrt(n1*p1*(1-p1)))
B <- matrix(c(p1,p12,p13,p14,p15,p12,p2,p23,p24,p25,p13,p23,p3,p34,p35,p14,p24,p34,p4,p45,p15,p25,p35,p45,p5),5,5,byrow = TRUE)
Theta1.1 <- matrix(rep(c(1,0,0,0,0),n1*alpha1),n1*alpha1,5,byrow = TRUE)
Theta1.2 <- matrix(rep(c(0,1,0,0,0),n1*alpha2),n1*alpha2,5,byrow = TRUE)
Theta1.3 <- matrix(rep(c(0,0,1,0,0),n1*alpha3),n1*alpha3,5,byrow = TRUE)
Theta1.4 <- matrix(rep(c(0,0,0,1,0),n1*alpha4),n1*alpha4,5,byrow = TRUE)
Theta1.5 <- matrix(rep(c(0,0,0,0,1),n1*alpha5),n1*alpha5,5,byrow = TRUE)
Theta1 <- rbind(Theta1.1,Theta1.2,Theta1.3,Theta1.4,Theta1.5)
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
extrempoint <- 2*sqrt(max(alpha1*psi1,alpha2*psi2,alpha3*psi3,alpha4*psi4,alpha5*psi5)+max((1-alpha1)*psi0,(1-alpha2)*psi0,(1-alpha3)*psi0,(1-alpha4)*psi0,(1-alpha5)*psi0))
###############################################################################################################################
#histogram of eigenvalues of A_theota
hist(e1,xlim=c(-1.2,1.2),ylim = c(0,45), xlab= NULL, ylab = "Frequency", probability = FALSE,breaks = 30,main=NULL)# 
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)


#histogram of eigenvalues of A_hat
hist(e1,xlim=c(-1.5,4),ylim = c(0,45), xlab= NULL, ylab = "Frequency", probability = FALSE,breaks = 40,main=NULL)# 
#lines(x_var2,y_var2,add=TRUE,col ="blue",lwd=2)
points(e2[1],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[2],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[3],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[4],0,col="black", lwd=10,pch=16,add=TRUE)
points(e2[5],0,col="black", lwd=10,pch=16,add=TRUE)
points(-extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
points(extrempoint,0,col="red", lwd=20,pch=16,add=TRUE)
##################################################################################




















#Run 1000 simulations

A_hat <- function(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5){
  n1
  p1
  p2
  p3
  p4
  p5
  p12
  p13
  p14
  p15
  p23
  p24
  p25
  p34
  p35
  p45
  alpha1
  alpha2
  alpha3
  alpha4
  alpha5
  gamma <- (1/sqrt(n1*p1*(1-p1)))
  B <- matrix(c(p1,p12,p13,p14,p15,p12,p2,p23,p24,p25,p13,p23,p3,p34,p35,p14,p24,p34,p4,p45,p15,p25,p35,p45,p5),5,5,byrow = TRUE)
  Theta1.1 <- matrix(rep(c(1,0,0,0,0),n1*alpha1),n1*alpha1,5,byrow = TRUE)
  Theta1.2 <- matrix(rep(c(0,1,0,0,0),n1*alpha2),n1*alpha2,5,byrow = TRUE)
  Theta1.3 <- matrix(rep(c(0,0,1,0,0),n1*alpha3),n1*alpha3,5,byrow = TRUE)
  Theta1.4 <- matrix(rep(c(0,0,0,1,0),n1*alpha4),n1*alpha4,5,byrow = TRUE)
  Theta1.5 <- matrix(rep(c(0,0,0,0,1),n1*alpha5),n1*alpha5,5,byrow = TRUE)
  Theta1 <- rbind(Theta1.1,Theta1.2,Theta1.3,Theta1.4,Theta1.5)
  ##4: generate matrix P=Theta*B*Theta^T
  #create matrix P
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

#eigenvalue 5
list_eigenvalue5_SBM5 <- c()
for (i in 1:1000) {
  list_eigenvalue5_SBM5[i] <- eigen(A_hat(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5),symmetric = TRUE)$values[3]
}

list_eigenvalue5_SBM5

sum(list_eigenvalue5_SBM5 >extrempoint)

##########################
#eigenvalue 6
list_eigenvalue6_SBM5 <- c()
for (i in 1:1000) {
  list_eigenvalue6_SBM5[i] <- eigen(A_hat(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5),symmetric = TRUE)$values[4]
}

list_eigenvalue6_SBM5

#list_eigenvalue_SBM2[5]

sum(list_eigenvalue6_SBM5 >extrempoint)
write.csv(list_eigenvalue6_SBM5, "list_eigenvalue6_SBM5")
write.csv(list_eigenvalue5_SBM5, "list_eigenvalue5_SBM5")
