library(Matrix)
library(klaR)
library(sbm)
# Set parameters
set.seed(88)
A <- function(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5){
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
  B <- matrix(c(p1,p12,p13,p12,p2,p23,p13,p23,p3),3,3,byrow = TRUE)
  ##3: Create matrix Theta
  Theta1.1 <- matrix(rep(c(1,0,0),n1*alpha1),n1*alpha1,3,byrow = TRUE)
  Theta1.2 <- matrix(rep(c(0,1,0),n1*alpha2),n1*alpha2,3,byrow = TRUE)
  Theta1.3 <- matrix(rep(c(0,0,1),n1*alpha3),n1*alpha3,3,byrow = TRUE)
  dim(Theta1.1)
  Theta1 <- rbind(Theta1.1,Theta1.2,Theta1.3)
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
  return(A)
}

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

K <- 5
n <- c(200, 200, 200, 200, 200)
block_sizes <- rep(n, each = K)
# Normalize adjacency matrix
hat_A <- A_hat(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5)
# Compute largest eigenvalues and eigenvectors
eigen_values <- eigen(hat_A, symmetric = TRUE)$values
eigen_vectors <- eigen(hat_A, symmetric = TRUE)$vectors

plot(eigen_values)
# Cluster the network using k-means
K_largest <- 5
embedding_vectors <- eigen_vectors[, 1:K_largest]
embedding_vectors <- apply(embedding_vectors, 2, function(x) x/sqrt(sum(x^2)))
X <- embedding_vectors
kmeans_clusters <- kmeans(X, K_largest)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score <- NMI(true_clusters, kmeans_clusters$cluster)

# Print normalized mutual information score
print(paste0("Normalized Mutual Information score: ", nmi_score))
colors <- c("red", "blue", "darkgreen", "gold1", "deeppink")
plot(X[,4:5],col = colors[true_clusters], pch = 16,xlab = " ",ylab = " ", main = "")
##########################



#plot(X_2,color=kmeans_clusters_2$cluster)






#Plot of SBM networks
library(igraph)
graph <- graph_from_adjacency_matrix(
  A(p1,p2,p3,p4,p5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,alpha1,alpha2,alpha3,alpha4,alpha5),
  mode = c( "undirected"),
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
d <- degree.distribution(graph)


plot(graph,main=("SBM graph with 5 communities"))


###############################################
#plot 3 eigenvectors
require(plotly)

data <- cbind(X,n)

data1 <- data.frame(data)


p <- plot_ly(data1, x=~X[,2], y=~X[,3], 
             z=~X[,4], color=~true_clusters,,size=10) %>%
  #add_markers(size=6) 
  layout(
    scene = list(
      xaxis = list(title = "Eigenvector 2"),
      yaxis = list(title = "Eigenvector 3"),
      zaxis = list(title = "Eigenvector 4")
    )
  )  
p <- p %>% hide_colorbar()
print(p)


