library(Matrix)
library(klaR)
library(sbm)
# Set parameters
A <- function(n1,p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3){
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

A_hat <- function(n1,p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3){
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

K <- 3
n <- c(300, 500, 200)
block_sizes <- rep(n, each = K)
# Normalize adjacency matrix
hat_A <- A_hat(n1,p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3)
# Compute largest eigenvalues and eigenvectors
eigen_values <- eigen(hat_A, symmetric = TRUE)$values
eigen_vectors <- eigen(hat_A, symmetric = TRUE)$vectors

plot(eigen_values)
# Cluster the network using k-means
K_largest <- 3
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
plot(X,color=true_clusters)
##########################
#Plot of SBM networks
library(igraph)
graph <- graph_from_adjacency_matrix(
  A(n1,p1,p2,p3,p12,p13,p23,alpha1,alpha2,alpha3),
  mode = c( "undirected"),
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
d <- degree.distribution(graph)


plot(graph,main=("SBM graph with 3 communities"))


###############################################
#plot 3 eigenvectors
require(plotly)

data <- cbind(X,n)

data1 <- data.frame(data)


p <- plot_ly(data1, x=~X[,1], y=~X[,2], 
             z=~X[,3], color=~true_clusters) %>%
  add_markers(size=6) 
layout(
    scene = list(
      xaxis = list(title = "Eigenvector 1"),
      yaxis = list(title = "Eigenvector 2"),
      zaxis = list(title = "Eigenvector 3")
  )
)  
p <- p %>% hide_colorbar()
print(p)


colors <- c("red", "blue", "darkgreen")
plot(X[,1:2],col = colors[true_clusters], pch = 16,xlab = " ",ylab = " ", main = "")
