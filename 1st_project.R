################################################################# 
#####            Project1: Bayesian Networks                #####
#################################################################
#
# 0) Setup environment
# 
### Install the pakages "deal" and "yeastCC"
install.packages("deal",repos="http://lib.stat.cmu.edu/R/CRAN",dependencies=TRUE)
BiocManager::install("yeastCC")

### Load the packages
library("deal")
library("yeastCC")
library("VGAM")

#
# 1) Data Preprocessing
#
### Load the data
expr <- as.data.frame(t(exprs(yeastCC)[orf800,]))
cat("Observations:", nrow(expr), "\n")
cat("Genes:", ncol(expr), "\n")




### Replace missing data with gene expression median
# for-loop
for (j in 1:ncol(expr)) {
  index_na <- which(is.na(expr[,j]))  # indices of NAs
  expr[index_na, j] <- median(expr[,j], na.rm=T)  # replace by median
}

# apply() alternative 
medianfill <- function(exprCol) {
  exprCol[which(is.na(exprCol))] <- median(exprCol, na.rm=T)
  return(exprCol)
}
expr <- as.data.frame(apply(expr, 2, medianfill))

### Filter genes based on inter quartile range iqr
# for-loop
iqr <- vector(length=ncol(expr))
for (j in 1:ncol(expr)) {
  iqr[j] <- quantile(expr[,j], c(0.75)) - quantile(expr[,j], c(0.25))
}

# apply() alternative
gene.iqr <- function(exprCol) {
  quantile(exprCol, c(0.75)) - quantile(exprCol, c(0.25))
}
iqr <- apply(expr, 2, gene.iqr)

# Keep only genes with iqr > 1.6
expr = expr[, iqr > 1.6]   # keep only genes with large variation
names(expr)  # print selected genes
#YBR054W	YRO2
#YBR088C	POL30
#YER124C	DSE1
#YGL028C	SCW11
#YLR286C	CTS1
#YHR143W	DSE2
#YNL327W	EGT2
#YGR108W	CLB1
#YNR067C	DSE4
#YOL007C	CSI2

# deal package seems not to work if all nodes are continuous.
# Thus, we need to add a discrete "dummy" node to work around this issue. 
# Adding this node will not have any influence on the final results.
expr$dummy <- factor(rep("1", nrow(expr)))
head(expr)

#
# 2-6) Build Bayesian Network
#
### 2) Create prior structure 
# EITHER: Using default 
#G0 <- network(expr) # prior structure with no edges

# OR: Specify your own prior network manually:
# To insert an arrow from node 'A' to node 'B',
# first click node 'A' and then click node 'B'.
# Add an arrow from YOL007C to YBR088C, 
# and from YNL327W to YER124C, YHR143W and YNR067C.
# When the graph is finished, click 'stop',
# Then, inspect the local probability distribution of node i,
# by clicking on node i and then on the background.,
# Finish by clicking the center button 'Stop'.
G0  <- network(expr, specifygraph=TRUE, inspectprob=TRUE)


# We don't want any arrows starting from the "dummy" node, thus we construct a list of banned dependencies:
banlist(G0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
plot(G0)


### 3) Show local probability distribution
localprob(G0)
localprob(G0)$YBR088C

### 4) Compute joint prior distribution
prior0 <- jointprior(G0, 5)  # equivalent to imaginary sample size = 5

### 5) Learn the initial network

G1 <- getnetwork(learn(G0, expr, prior0))
print(G1$score)

### 6) Search for optimal network (takes some time)
nwSearch <- autosearch(G1, expr, prior0, removecycles=FALSE, trace=FALSE)
G <- getnetwork(nwSearch)
plot(G)

### Function for building an optimal network from expression data (from above)
build.optimal.network <- function(exprData) { 
  N0 <- getnetwork(learn(G0, exprData, prior0))
  nwasarch <- autosearch(N0, exprData, prior0, removecycles=FALSE, trace=FALSE)
  getnetwork(nwasarch)
}

### Custom function for plotting a BN 
plot.bn <- function(BN, file=NULL) {
  par(mar=c(0,0,0,0))
  plot(BN, cexscale=13, unitscale=27, arrowlength=0.1, xr=c(0, 350), yr=c(20,370))
  if (!is.null(file)) {
    plt <- recordPlot()	
    pdf(file)
    replayPlot(plt)
    dev.off()
  }
}
BN <- build.optimal.network(expr)
plot(BN)
plot.bn(BN, file="p_BNstar5.pdf")

genes = within(expr, rm(dummy))
genes.vars = sapply(genes, var)

perturbed_data = list()
for (i in 1:30) {
  p_genes = data.frame(genes)
  for (gene in colnames(p_genes)) {
    p_genes[gene] = p_genes[gene] + rnorm(nrow(p_genes), mean=0, sd=genes.vars[gene] / 10)
  }
  perturbed_data[[i]] = p_genes
}

library(reshape2)
library(ggplo2)

yhr143w = t(sapply(perturbed_data, function(x) x$YHR143W))
yhr143w = data.frame(yhr143w)
colnames(yhr143w) = 1:77
yhr143w.melted = melt(yhr143w)
ggplot(yhr143w.melted, aes(x=variable, y=value)) + geom_boxplot()

for (i in 1:30) {
  perturbed_data[[i]]$dummy = factor(rep("1", nrow(perturbed_data[[i]])))
}

p_networks = list()
for (i in 1:30) {
  p_networks[[i]] = build.optimal.network(perturbed_data[[i]])
}
plot(p_networks[[5]])
plot.bn(p_networks[[5]], file="PBN5.pdf")

get_edges = function(network) {
  edges = list()
  i = 1
  for (node in network$nodes) {
    for (parent_node in node$parents) {
      edges[[i]] = c(parent_node, node$idx)
      i = i + 1
    }
  }
  return(edges)
}

df_edges_freqs = function(net, edges, all_edges) {
  edges_freqs = list()
  for (edge in edges) {
    count = 0
    for (p_edge in all_edges) {
      if (identical(edge, p_edge)) {
        count = count + 1
      }
    }
    freq = count / 30
    edges_freqs = c(edges_freqs, freq)
  }
  
  edges_names = list()
  for (edge in edges) {
    edge_name = paste(net$nodes[[edge[1]]]$name, net$nodes[[edge[2]]]$name, sep="->")
    edges_names = c(edges_names, edge_name)
  }
  
  edges_df = data.frame(freq=unlist(edges_freqs), name=unlist(edges_names))
  
  return(edges_df)
}

plot_edges_df = function(edges_df) {
  g = ggplot(edges_df, aes(x=name, y=freq)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(g)
}

BN_edges = get_edges(BN)
p_networks_edges = list()
for (i in 1:30) {
  p_networks_edges = c(p_networks_edges, get_edges(p_networks[[i]]))
}

BN_edges_df = df_edges_freqs(BN, BN_edges, p_networks_edges)
BN_edges_df[BN_edges_df$freq < 0.5,]
plot_edges_df(BN_edges_df)

unique_edges = p_networks_edges[!duplicated(p_networks_edges)]
unique_edges_wo_BN = setdiff(unique_edges, BN_edges)

not_BN_edges_df = plot_edges_freqs(BN, unique_edges_wo_BN, p_networks_edges)
not_BN_edges_df[not_BN_edges_df$freq >= 0.33,]
plot_edges_df(not_BN_edges_df)
