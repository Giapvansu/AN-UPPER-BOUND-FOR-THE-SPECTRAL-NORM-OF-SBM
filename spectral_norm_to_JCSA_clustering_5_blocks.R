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

################################################

################################################
################################################
####NMI for 2 clusters
# Cluster the network into two clusters using k-means
K_2 <- 2
embedding_vectors_2 <- eigen_vectors[, 1:K_2]
embedding_vectors_2 <- apply(embedding_vectors_2, 2, function(x) x/sqrt(sum(x^2)))
X_2 <- embedding_vectors_2
kmeans_clusters_2 <- kmeans(X_2, K_2)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_2 <- NMI(true_clusters, kmeans_clusters_2$cluster)

# Print normalized mutual information score two blocks
print(paste0("Normalized Mutual Information score two clusters: ", nmi_score_2))
colors2 <- c("red", "blue", "darkgreen")
plot(X_2,col = colors2[kmeans_clusters_2$cluster], pch = 16,xlab = " ",ylab = " ", main = "")

################################
#NMI for 3 clusters
# Cluster the network into three clusters using k-means
K_3 <- 3
embedding_vectors_3<- eigen_vectors[, 1:K_3]
embedding_vectors_3 <- apply(embedding_vectors_3, 3, function(x) x/sqrt(sum(x^2)))
X_3 <- embedding_vectors_3
kmeans_clusters_3 <- kmeans(X_3, K_3)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_3 <- NMI(true_clusters, kmeans_clusters_3$cluster)

####NMI for 4 clusters
# Cluster the network into 4 clusters using k-means
K_4 <- 4
embedding_vectors_4 <- eigen_vectors[, 1:K_4]
embedding_vectors_4 <- apply(embedding_vectors_4, 2, function(x) x/sqrt(sum(x^2)))
X_4 <- embedding_vectors_4
kmeans_clusters_4 <- kmeans(X_4, K_4)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_4 <- NMI(true_clusters, kmeans_clusters_4$cluster)

# Print normalized mutual information score 4 blocks
print(paste0("Normalized Mutual Information score four clusters: ", nmi_score_4))
colors4 <- c("red", "blue", "darkgreen","orange")
plot(X_4[,1:2],col = colors4[kmeans_clusters_4$cluster], pch = 16,xlab = " ",ylab = " ", main = "")


######################
################################################
####NMI for 5 clusters
# Cluster the network into 5 clusters using k-means
K_5 <- 5
embedding_vectors_5 <- eigen_vectors[, 1:K_5]
embedding_vectors_5 <- apply(embedding_vectors_5, 2, function(x) x/sqrt(sum(x^2)))
X_5 <- embedding_vectors_5
kmeans_clusters_5 <- kmeans(X_5, K_5)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_5 <- NMI(true_clusters, kmeans_clusters_5$cluster)

# Print normalized mutual information score 5 blocks
print(paste0("Normalized Mutual Information score five clusters: ", nmi_score_5))
colors5 <- c("red", "blue", "darkgreen","orange","aquamarine")
plot(X_5[,1:2],col = colors5[kmeans_clusters_5$cluster], pch = 16,xlab = " ",ylab = " ", main = "")


#################################################
####NMI for 6 clusters
# Cluster the network into 5 clusters using k-means
K_6 <- 6
embedding_vectors_6 <- eigen_vectors[, 1:K_6]
embedding_vectors_6 <- apply(embedding_vectors_6, 2, function(x) x/sqrt(sum(x^2)))
X_6 <- embedding_vectors_6
kmeans_clusters_6 <- kmeans(X_6, K_6)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_6 <- NMI(true_clusters, kmeans_clusters_6$cluster)

# Print normalized mutual information score 6 blocks
print(paste0("Normalized Mutual Information score six clusters: ", nmi_score_6))
colors6 <- c("red", "blue", "darkgreen","orange","aquamarine","black")
plot(X_6[,1:2],col = colors6[kmeans_clusters_6$cluster], pch = 16,xlab = " ",ylab = " ", main = "")


#################################################
####NMI for 7 clusters
# Cluster the network into 7 clusters using k-means
K_7 <- 7
embedding_vectors_7 <- eigen_vectors[, 1:K_7]
embedding_vectors_7 <- apply(embedding_vectors_7, 2, function(x) x/sqrt(sum(x^2)))
X_7 <- embedding_vectors_7
kmeans_clusters_7 <- kmeans(X_7, K_7)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_7 <- NMI(true_clusters, kmeans_clusters_7$cluster)

# Print normalized mutual information score 6 blocks
print(paste0("Normalized Mutual Information score seven clusters: ", nmi_score_7))
colors7 <- c("red", "blue", "darkgreen","orange","aquamarine","black", "pink")
plot(X_7[,1:2],col = colors7[kmeans_clusters_7$cluster], pch = 16,xlab = " ",ylab = " ", main = "")


####NMI for 8 clusters
# Cluster the network into 8 clusters using k-means
K_8 <- 8
embedding_vectors_8 <- eigen_vectors[, 1:K_8]
embedding_vectors_8 <- apply(embedding_vectors_8, 2, function(x) x/sqrt(sum(x^2)))
X_8 <- embedding_vectors_8
kmeans_clusters_8 <- kmeans(X_8, K_8)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_8 <- NMI(true_clusters, kmeans_clusters_8$cluster)

# Print normalized mutual information score 8 blocks
print(paste0("Normalized Mutual Information score seven clusters: ", nmi_score_8))
colors8 <- c("red", "blue", "darkgreen","orange","aquamarine","black", "pink","darkred")
plot(X_8[,1:2],col = colors8[kmeans_clusters_8$cluster], pch = 16,xlab = " ",ylab = " ", main = "")


####NMI for 9 clusters
# Cluster the network into 9 clusters using k-means
K_9 <- 9
embedding_vectors_9 <- eigen_vectors[, 1:K_9]
embedding_vectors_9 <- apply(embedding_vectors_9, 2, function(x) x/sqrt(sum(x^2)))
X_9 <- embedding_vectors_9
kmeans_clusters_9 <- kmeans(X_9, K_9)
# Compute normalized mutual information score
library(aricode)
true_clusters <- rep(1:K, times = n)
nmi_score_9 <- NMI(true_clusters, kmeans_clusters_9$cluster)

# Print normalized mutual information score 8 blocks
print(paste0("Normalized Mutual Information score seven clusters: ", nmi_score_9))
colors9 <- c("red", "blue", "darkgreen","orange","aquamarine","black", "pink","darkred","green")
plot(X_9[,1:2],col = colors9[kmeans_clusters_9$cluster], pch = 16,xlab = " ",ylab = " ", main = "")





hist(rowSums(hat_A)/n1)

#plot of NMI score
K_x <- c(2,3,4,5,6,7,8,9)
nmi_list <- c(nmi_score_2,nmi_score_3,nmi_score_4,nmi_score,nmi_score_6,nmi_score_7,nmi_score_8,nmi_score_9)
plot(K_x,nmi_list,type = "b", pch = 19, 
     col = "red", xlab="number of clusters", ylab="NMI score",main="NMI score")


