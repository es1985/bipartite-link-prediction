require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)
require(SnowballC)
require(lsa)
require(scales)
require(reshape2)
options(java.parameters = "-Xmx4g")

rm(list=ls(all=TRUE))

#CAIDA
edgelist_source <- read.csv('/home/esavin/Link prediction/as-caida20071029.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))
edgelist_target <- read.csv('/home/esavin/Link prediction/as-caida20071105.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))
edgelist_pred <- read.csv('/home/esavin/Link prediction/as-caida20071112.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))

#UCLA
edgelist_source <- read.csv('/home/esavin/Link prediction/UCLA/201310.relationship',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))
edgelist_target <- read.csv('/home/esavin/Link prediction/UCLA/201403.relationship',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))
edgelist_pred <- read.csv('/home/esavin/Link prediction/UCLA/201410.relationship',sep = "\t",header = TRUE, row.names = NULL, col.names = c('as1','as2','relationship'))

edgelist_target <- unique(rbind(edgelist_source,edgelist_target))
edgelist_pred <-unique(rbind(edgelist_target,edgelist_pred))

#Equate the number of nodes in source, target and test datasets
outersect <- function (x,y)
{
  sort(c(setdiff(x,y),setdiff(y,x)))
}

diff_as1<-outersect(edgelist_source$as1,edgelist_target$as1)
diff_as2<-outersect(edgelist_source$as2,edgelist_target$as2)

while(length(diff_as1) > 0 & length(diff_as2) > 0)
{
  for (i in 1:length(diff_as1))
  {
    edgelist_source <- edgelist_source[edgelist_source$as1 != diff_as1[i],]
    edgelist_target <- edgelist_target[edgelist_target$as1 != diff_as1[i],]
  }
  
  for (i in 1:length(diff_as2))
  {
    edgelist_source <- edgelist_source[edgelist_source$as2 != diff_as2[i],]
    edgelist_target <- edgelist_target[edgelist_target$as2 != diff_as2[i],]
  }
  
  
  diff_as1<-outersect(edgelist_source$as1,edgelist_pred$as1)
  diff_as2<-outersect(edgelist_source$as2,edgelist_pred$as2)
  
  #Remove rows where elements appear in just 1 data set and not both, i.e. make sure both datasets
  #have the same number of nodes
  for (i in 1:length(diff_as1))
  {
    edgelist_source <- edgelist_source[edgelist_source$as1 != diff_as1[i],]
    edgelist_pred <- edgelist_pred[edgelist_pred$as1 != diff_as1[i],]
  }
  
  for (i in 1:length(diff_as2))
  {
    edgelist_source <- edgelist_source[edgelist_source$as2 != diff_as2[i],]
    edgelist_pred <- edgelist_pred[edgelist_pred$as2 != diff_as2[i],]
  }
  
  diff_as1<-outersect(edgelist_source$as1,edgelist_target$as1)
  diff_as2<-outersect(edgelist_source$as2,edgelist_target$as2)
}


#Detect new edges between the target and the source data sets
#new_edges_ts <-edgelist_target[!duplicated(rbind(edgelist_source, edgelist_target))[-seq_len(nrow(edgelist_source))], ]

#Filter out as1s active in the source dataset 
#active_as1s_source <-edgelist_source[edgelist_source$as1 %in% new_edges_ts$as1,]

#Detect new edges between the prediction and target
#new_edges_pt <-edgelist_pred[!duplicated(rbind(edgelist_target, edgelist_pred))[-seq_len(nrow(edgelist_target))], ]

#Filter out as1s active in the source dataset
#active_as1s_target <-edgelist_target[edgelist_target$as1 %in% new_edges_pt$as1,]

sample_a <- edgelist_source
sample_b <- edgelist_target
sample_c <- edgelist_pred


#set the number of eigenvalues to be caluclated

r=100

g_a <- graph.data.frame(sample_a,directed = FALSE)
#is.connected(g_a)
#cl<-clusters(g_a)
#g_a <- subgraph(g_a, which(cl$membership == which.max(cl$csize)))
#is.connected(g_a)
#vcount(g_a)
adj_a <- get.adjacency(g_a, sparse = TRUE,type=c("both"))
#image(adj_a)

#calculate the eigenvalues and eigenvectors using ARPACK algorithm



g_b <- graph.data.frame(sample_b,directed = FALSE)
#is.connected(g_b)
#cl<-clusters(g_b)
#g_b <- subgraph(g_b, which(cl$membership == which.max(cl$csize)))
#is.connected(g_b)
#vcount(g_b)
adj_b <- get.adjacency(g_b, sparse = TRUE)
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_b <- adj_b[row.order,col.order]
#image(adj_b)

#sanity check - reconstruct the original matrix A
#lamda = Matrix(0,  ncol = r,nrow = r)
#diag(lamda) = baev_a$values
#t<- baev_a$vectors %*% lamda %*% t(baev_a$vectors)
#hist(as.numeric(t))
#t<- rescale(t,to = c(0,1))
#hist(as.numeric(t))
#substitute values in a matrix
#t[t<0.3] <- 0
#t[t>0.3] <- 1
#t<-Matrix(t,sparse = TRUE)
#image(adj_a)
#image(t)


#Check the development of the prediction eigenvalues in comparison to target
g_c <- graph.data.frame(sample_c,directed = FALSE)
#V(g_c)$type <- V(g_c)$name %in% sample_c[,1]
#get.edgelist(g_c) #to verify this was done properly
adj_c <- get.adjacency(g_c, sparse = TRUE)
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_c <- adj_c[row.order,col.order]


f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a %*% x) }
baev_a <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b %*% x) }
baev_b <- arpack(f2, sym=TRUE, options=list(n=vcount(g_b), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_b)*12))

#f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_c %*% x) }
#baev_c <- arpack(f2, sym=TRUE, options=list(n=vcount(g_c), nev=r, ncv=r+3,
#                                            which="LM", maxiter=vcount(g_c)*15))


#compute matrix of new edges - i.e. B
ne <- adj_b - adj_a
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
ne <- ne[row.order,col.order]

#image(adj_a)
#image(adj_b)
#image(ne)


#normalise the spectra
nlamda <- baev_a$values * abs(baev_a$values[1])/abs(baev_b$values[1])

#without normalisation
#nlamda <- baev_a$values


#compute multiplication of eigenvalues of A and matrix B
b_eigen <- t(baev_a$vectors) %*% ne %*% baev_a$vectors
#hist(as.numeric(b_eigen))
image(b_eigen)


#df = as.data.frame(nlamda)
#df = as.data.frame(as.numeric(diag(b_eigen))) #the independent variable
#df$flamda = nlamda #the dependent function
#colnames(df) = c("source","target")
#plot(df$source,df$target)

df = as.data.frame(baev_a$values) #the independent variable
df$flamda = as.numeric(diag(b_eigen)) #the dependent function
colnames(df) = c("source","target")
plot(x=df$source,y=df$target)

#path counting
a_polyn <-nls(formula = target ~ I(alpha*source) + I(beta*source^2) + I(gamma*source^3) + I(delta*source^4) + I(delta*source^5),start = c(alpha=-10,beta=-10,gamma=-10,delta=-10), data = df)
summary(a_polyn)
sum(resid(a_polyn)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_polyn,new)
points(x,y,col = "purple")

#nonnegative polynomial path counting
a_polyn <-nls(formula = target ~ I(alpha*source) + I(beta*source^2) + I(gamma*source^3) + I(delta*source^4) + I(delta*source^5),start = c(alpha=0.01,beta=0.01,gamma=0.01,delta=0.01), data = df)
summary(a_polyn)
sum(resid(a_polyn)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_polyn,new)
points(x,y,col = "green")

#matrix exponential
a_exp <- nls(formula = target ~ I(alpha*source)+I(((alpha*source)^2)/2)+I(((alpha*source)^3)/6+I(((alpha*source)^4)/24)),start=c(alpha=0.01),data=df)
summary(a_exp)
sum(resid(a_exp)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_exp,new)
points(x,y,col = "red")

#Neumann pseudokernel
a_neu <- nls(formula = target ~ I(alpha*source)+I((alpha*source)^2)+I((alpha*source)^3)+I((alpha*source)^4)+I((alpha*source)^5),start=c(alpha=0.01),data = df)
summary(a_neu)
sum(resid(a_neu)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_neu,new)
points(x,y,col = "blue")


#Now construct a prediction based on the best fit
#Now construct a prediction based on the best fit
flamda_a_polyn = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_polyn) = predict(a_polyn,new.dfb)
mp_a_polyn <- baev_b$vectors %*% flamda_a_polyn %*% t(baev_b$vectors)
dimnames(mp_a_polyn) <- dimnames(adj_a)

flamda_a_exp = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_exp) = predict(a_exp,new.dfb)
mp_a_exp <- baev_b$vectors %*% flamda_a_exp %*% t(baev_b$vectors)
dimnames(mp_a_exp) <- dimnames(adj_a)

flamda_a_neu = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_neu) = predict(a_neu,new.dfb)
mp_a_neu <- baev_b$vectors %*% flamda_a_neu %*% t(baev_b$vectors)
dimnames(mp_a_neu) <- dimnames(adj_a)



#
netest <- adj_c - adj_b
hist(as.numeric(netest))
image(netest) #new edges
#substitute any multiple edges with just 1 edge
netest[netest >= 1] <- 1
netest[netest < 0] <- 0

pred_a_polyn <- prediction(as.vector(mp_a_polyn),as.vector(netest))

perf_a_polyn <- performance(pred_a_polyn, "tpr", "fpr")
plot(perf_a_polyn)
precision_recall_a_polyn <- performance(pred_a_polyn, "prec", "rec")
plot(precision_recall_a_polyn)
sensitivity_specificity_a_polyn <- performance(pred_a_polyn,"sens","spec")
plot(sensitivity_specificity_a_polyn)
lift_a_polyn <- performance(pred,"lift","rpp")
plot(lift_a_polyn)

acc_pred_a_polyn <- performance(pred_a_polyn,'acc')
f_a_polyn <- performance(pred_a_polyn,'f')
plot(f_a_polyn)
auc_a_polyn <- performance(pred_a_polyn,"auc")


pred_a_exp <- prediction(as.vector(mp_a_exp),as.vector(netest))

perf_a_exp <- performance(pred_a_exp, "tpr", "fpr")
plot(perf_a_exp)
precision_recall_a_exp <- performance(pred_a_exp, "prec", "rec")
plot(precision_recall_a_exp)
sensitivity_specificity_a_exp <- performance(pred_a_exp,"sens","spec")
plot(sensitivity_specificity_a_exp)
lift_a_exp <- performance(pred,"lift","rpp")
plot(lift_a_exp)

acc_pred_a_exp <- performance(pred_a_exp,'acc')
err_pred_a_exp <- performance(pred_a_exp,'err')
f_a_exp <- performance(pred_a_exp,'f')
plot(f_a_exp)
auc_a_exp <- performance(pred_a_exp,"auc")


pred_a_neu <- prediction(as.vector(mp_a_neu),as.vector(netest))

perf_a_neu <- performance(pred_a_neu, "tpr", "fpr")
plot(perf_a_neu)
precision_recall_a_neu <- performance(pred_a_neu, "prec", "rec")
plot(precision_recall_a_neu)
sensitivity_specificity_a_neu <- performance(pred_a_neu,"sens","spec")
plot(sensitivity_specificity_a_neu)
lift_a_neu <- performance(pred,"lift","rpp")
plot(lift_a_neu)

acc_pred_a_neu <- performance(pred_a_neu,'acc')
f_a_neu <- performance(pred_a_neu,'f')
plot(f_a_neu)
auc_a_neu <- performance(pred_a_neu,"auc")

#calculate jaccard similarity measure
adj_b_jacc <- Matrix(similarity.jaccard(g_b),sparse=TRUE)

pred_jacc <- prediction(as.vector(adj_b_jacc),as.vector(netest))
perf_jacc <- performance(pred,"tpr","fpr")
plot(perf_jacc)
sensitivity_specificity_jacc <- performance(pred,"sens","spec")
plot(sensitivity_specificity_jacc)
precision_recall_jacc <- performance(pred, "prec", "rec")
plot(precision_recall_jacc)
auc_jacc <- performance(pred_jacc,"auc")
auc_jacc

#calculate adamic adar similarity measure
adj_b_adar <- Matrix(similarity.invlogweighted(g_b),sparse=TRUE)

pred_adar <- prediction(as.vector(adj_b_adar),as.vector(netest))
perf_adar <- performance(pred,"tpr","fpr")
plot(perf_adar)
sensitivity_specificity_adar <- performance(pred,"sens","spec")
plot(sensitivity_specificity_adar)
precision_recall_adar <- performance(pred_adar, "prec", "rec")
plot(precision_recall_adar)
auc_adar <- performance(pred_adar,"auc")
auc_adar


degree <- function(m)
  #The function accepts adjacency matrix and returns the degree matrix
{
  D <- Matrix(0,nrow = nrow(m),ncol = nrow(m))
  for (i in 1:nrow(m))
  { 
    # print(paste('D calculation : row', i))
    D[i,i] = sum(m[i,])
  }
  return (D)
}

adj_b_d<-degree(adj_b)

#calculate preferential attachment scores
adj_b_ppa<-as.matrix(as.vector(diag(adj_b_d))%*%t(as.vector(diag(adj_b_d))))

pred_ppa <- prediction(as.vector(adj_b_ppa),as.vector(netest))
perf_ppa <- performance(pred_ppa,"tpr","fpr")
plot(perf_ppa)
sensitivity_specificity_ppa <- performance(pred,"sens","spec")
plot(sensitivity_specificity_ppa)
precision_recall_ppa <- performance(pred_ppa, "prec", "rec")
plot(precision_recall_ppa)
auc_ppa <- performance(pred_ppa,"auc")
auc_ppa


cutoffs <- data.frame(cut=perf_ppa@alpha.values[[1]], fpr=perf_ppa@x.values[[1]], 
                      tpr=perf_ppa@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
head(cutoffs)
head(subset(cutoffs, fpr < 0.50))

cutoffs <- data.frame(cut=precision_recall@alpha.values[[1]], recall=precision_recall@x.values[[1]], 
                      precision=precision_recall@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$precision, decreasing=TRUE),]
head(cutoffs)
head(subset(cutoffs, fpr < 0.20))

#normalise the matrix

adj_a_d <- degree(adj_a)

#take squre root of the diagonal of the degree matrix
adj_a_ds <- Matrix(data = 0,nrow=nrow(adj_a),ncol=nrow(adj_a))

diag(adj_a_ds)<-1/sqrt(diag(adj_a_d))

adj_a_n <- adj_a_ds%*%adj_a%*%adj_a_ds

t <- (adj_a_n + t(adj_a_n))

ne_n <- adj_b + t(adj_b) - (adj_a + t(adj_a))


f2 <- function(x, extra=NULL) { cat("."); as.vector((adj_a_n + t(adj_a_n)) %*% x) }
baev_a_n <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))

#compute multiplication of eigenvalues of A and matrix B
b_eigen_n <- t(baev_a_n$vectors) %*% ne %*% baev_a_n$vectors
#hist(as.numeric(b_eigen))
image(b_eigen_n)

df = as.data.frame(baev_a_n$values) #the independent variable
df$flamda = as.numeric(diag(b_eigen)) #the dependent function
colnames(df) = c("source","target")
plot(x=df$source,y=df$target)

#polynomial
n_poly <-nls(formula = target ~ I(alpha*source) + I(beta*source^2) + I(gamma*source^3) + I(delta*source^4) + I(delta*source^5),start = c(alpha=0.01,beta=0.01,gamma=0.01,delta=0.01), data = df)
summary(n_poly)
sum(resid(n_poly)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_poly,new)
points(x,y,col = "green")

#matrix exponential
n_exp <- nls(formula = target ~ I(alpha*source)+I(((alpha*source)^2)/2)+I(((alpha*source)^3)/6+I(((alpha*source)^4)/24)),start=c(alpha=0.01),data=df)
summary(n_exp)
sum(resid(n_exp)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_exp,new)
points(x,y,col = "red")

#Neumann pseudokernel
n_neu <- nls(formula = target ~ I(alpha*source)+I((alpha*source)^2)+I((alpha*source)^3)+I((alpha*source)^4)+I((alpha*source)^5),start=c(alpha=0.01),data = df)
summary(n_neu)
sum(resid(n_neu)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_neu,new)
points(x,y,col = "blue")



#######################################################################

#Compare the prediction and the real matrix
ind <- which(mp > 2.305964e-05)
for (i in 1:length(ind))
{
  k <- arrayInd(ind[i],dim(mp))
  rownames(mp)[k[,1]]
  colnames(mp)[k[,2]]
  print(paste('Predicted positive score',mp[ind[i]],' for the the following link'))
  print(mapply(`[[`, dimnames(mp), k))
  print(paste('In reality the following value was observed in the adjacency matrix'))
  print(adj_c[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
  print(paste('New edges matrix had the following value'))
  print(netest[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
}

#' This function computes the average precision at k
#' between two sequences
#'
#' @param k max length of predicted sequence
#' @param actual ground truth set (vector)
#' @param predicted predicted sequence (vector)
#' @export 

apk <- function(k, actual, predicted)
{
  score <- 0.0
  cnt <- 0.0
  for (i in 1:min(k,length(predicted)))
  {
    if (predicted[i] %in% actual && !(predicted[i] %in% predicted[0:(i-1)]))
    {
      cnt <- cnt + 1
      score <- score + cnt/i 
    }
  }
  score <- score / min(length(actual), k)
  score
}

#' Compute the mean average precision at k
#'
#' This function computes the mean average precision at k
#' of two lists of sequences.
#'
#' @param k max length of predicted sequence
#' @param actual list of ground truth sets (vectors)
#' @param predicted list of predicted sequences (vectors)
#' @export

mapk <- function (k, actual, predicted)
{
  scores <- rep(0, length(actual))
  for (i in 1:length(scores))
  {
    scores[i] <- apk(k, actual[[i]], predicted[[i]])
  }
  score <- mean(scores)
  score
}

apk(450,as.vector(netest),as.vector(mp))

apk(450,as.vector(netest),as.vector(mp))



#matrix normalisation
normalise <- function(m)
{
  D = Matrix(data = 0,nrow=nrow(m),ncol=nrow(m))
  #The function normalises the matrix as described on pages 15-16 
  for (i in 1:nrow(m))
  { 
    print(paste('D1 calculation : row', i))
    D1[i,i] = sum(m[i,])
  }
  D2 = Matrix(data = 0,nrow = ncol(m),ncol=ncol(m))  
  for (i in 1:ncol(m))
  {  
    print(paste('D2 calcilation : column',i))
    D2[i,i] = sum(m[,i])
  }
  diag(D1) <- 1/sqrt(diag(D1))
  diag(D2) <- 1/sqrt(diag(D2))
  M=D1%*%m%*%D2
  return (M)
}

adj_an <- normalise(adj_a)
hist(as.numeric(adj_an))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_an %*% x) }
baev_an <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_bn %*% x) }
baev_bn <- arpack(f2, sym=TRUE, options=list(n=vcount(g_b), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_b)*12))

bla <- t(adj_an)

df = as.data.frame(baev_an$values) #the independent variable
df$flamda = as.numeric(diag(b_eigen)) #the dependent function
colnames(df) = c("source","target")
plot(x=df$source,y=df$target)

#convert matrix of prediction into 0/1 matrix
h<-hist(as.numeric(mp))
#find concetration of values and set them to be 0
#z<-match(max(h$counts),h$counts) #zero 
h$counts
h$breaks
l<-length(h$counts)
#h$breaks[z]
#mp[mp<=h$breaks[1]] <- 0
mp[mp>=h$breaks[l-round(l/1.3)]] <- 1
mp[mp<h$breaks[l-round(l/1.3)]] <- 0
hist(as.numeric(mp))
image(mp)
netest <- adj_c - adj_b
hist(as.numeric(netest))
image(netest)


#Create confusion matrix
mx<-max(as.numeric(netest))
mn<-min(as.numeric(netest))

m<-netest+1.5*mp
hist(as.numeric(m))
dimnames(m) <- dimnames(adj_a)
#m<-2*as.matrix(mp)+as.matrix(netest)
##image(m) #new edges
#m<-2*as.matrix(mp)+as.matrix(netest)
#image(m) #new edges
tp <- which(m >= mx+1.5, arr.in = TRUE, useNames = TRUE)
fp <- which(m == 1.5, arr.in = TRUE, useNames = TRUE)
fn <- which(m == mx, arr.in = TRUE, useNames = TRUE)
tn <- which(m <= 0, arr.in = TRUE, useNames = TRUE)


confusion_matrix <- matrix(c(nrow(tp),nrow(fp),nrow(fn),nrow(tn)),nrow = 2,ncol = 2,byrow = TRUE)
confusion_matrix

precision = nrow(tp)/(nrow(tp)+nrow(fp)) #positive predictive value
npv = nrow(tn)/(nrow(fn)+nrow(tn)) #negative predictive value
recall = nrow(tp)/(nrow(fn)+nrow(tp)) #sensitivity
specificity = nrow(tn)/(nrow(fp)+nrow(tn))

precision_matrix <- matrix(c(precision,npv,recall,specificity),nrow = 2,ncol = 2,byrow = TRUE)
precision_matrix



##assemble back into the data frame
arrayInd(t,useNames=TRUE)

t
m[1790,5]

pred <- melt(m)
pred <- pred[,1:2]

pred <- pred[pred$value == 1,]
pred <- pred[order(pred[1,])]
preddf <- NULL
preddf <-as.data.frame(pred$Var2)
preddf$as2 <- pred$Var1
colnames(preddf) <- c("as1","as2")

dimnames(m==2)

m>0

m<-Matrix(m,sparse = TRUE)
image(m)

pred = as.matrix(as.numeric(m))


precision <- subset(m,subset = m[]==0)

require(pROC)
roc

pred <- as.matrix(m)
pred <- melt(pred)
pred <- pred[pred$value == 1,]
pred <- pred[,1:2]
pred <- pred[order(pred[1,])]
preddf <- NULL
preddf <-as.data.frame(pred$Var2)
preddf$as2 <- pred$Var1
colnames(preddf) <- c("as1","as2")

match <-new_edges_test[!duplicated(rbind(new_edges_test, preddf))[-seq_len(nrow(new_edges_test))], ]









###############################################################################################


#function for selecting a subsample of transacting customers and as2s from 3 different datasets
sample_a <- edgelist_source
sample_a <- NULL
sample_b <- edgelist_source
sample_b <- NULL
sample_c <- edgelist_source
sample_c <- NULL

x <- as.integer(runif(10,1,nrow(edgelist_source))) 
#x <- active_as1s_source[1:5,]
#d<-edgelist_source[1:5,]$as1
d<-edgelist_source[x,]$as1
for (i in 1 : length(d))
{
  #i=1
  sc_source<-edgelist_source[edgelist_source$as1==d[i],] #set of as1s source
  sc_target<-edgelist_target[edgelist_target$as1==d[i],] #set of as1s target
  sc_pred<-edgelist_pred[edgelist_pred$as1==d[i],] #set of as1s prediction
  
  cst<-setdiff(sc_target$as2,sc_source$as2) #check that there were no new as2s between source and target 
  ctp<-setdiff(sc_pred$as2,sc_source$as2) #check that there were no new as2s added between target and prediction 
  
  for (j in 1 : length(cst))
  {
    sc_target <- sc_target[sc_target$as2 != cst[j],]
  }
  
  for (l in 1 : length(ctp))
  {
    sc_pred <- sc_pred[sc_pred$as2 != ctp[l],]
  }
  
  #add all source rows to sample_a
  sample_a <- rbind(sample_a,sc_source)
  
  #add all source rows to sample_b also add new rows to sample_b
  sample_b <- rbind(sample_b,sc_target)
  
  #add all source rows, target rows and prediction rows to sample_c
  sample_c <- rbind(sample_c,sc_pred)
}



t<-sample_a[!duplicated(sample_a),] #remove duplicated edges

sample_b <- rbind(sample_b,sample_a)
t1<-sample_b[!duplicated(sample_b),] #remove duplicated edges

sample_c <- rbind(sample_c,sample_b)
t2 <- sample_c[!duplicated(sample_c),] #remove duplicated edges


#verify the as2 sets are identical
setdiff(sample_b$as2,sample_a$as2)
setdiff(sample_c$as2,sample_b$as2)

setdiff(sample_a$as2,sample_b$as2)
setdiff(sample_b$as2,sample_c$as2)

#sample_a<-sc_source #after data frame
#sample_b<-sc_target #before data frame
#sample_c<-sc_pred #predict data frame


fit1 <- nls(formula = df$target ~ sinh(alpha*df$source),start=c(alpha=0.01),data = df)
summary(fit1)
sum(resid(fit1)^2)
confint(fit1)
new = data.frame(xdata = seq(min(df$source),max(df$source),len=length(df$source)))
lines(df$source,predict(fit1,newdata=new),col = "red")


fit2 <- nls(formula = df$target ~ alpha*(I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9)), start=c(alpha=0.01),data = df)
summary(fit2)
sum(resid(fit2)^2)
confint(fit2)
new = data.frame(xdata = seq(min(df$source),max(df$source),len=length(df$source)))
#lines(x=new$xdata,y=predict(fit2,newdata=new))
lines(df$source,y=predict(fit2,newdata=new),col = "green")


fit3 <- nls(formula = df$target ~ alpha*df$source/(1-(alpha*df$source)^2), start=c(alpha=0.1),data = df)
summary(fit3)
sum(resid(fit3)^2)
confint(fit3)
new = data.frame(xdata = seq(min(df$source),max(df$source),len=length(df$source)))
lines(df$source,predict(fit3,newdata=new),col="blue")

#alpha <- 8.714e-11
#df$source <- alpha*df$source
#ggplot(data=df,aes(source,target))+geom_point()+
#  stat_smooth(method = "lm", formula = df$target ~ (I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9)),col = 'blue')+
#  stat_smooth(method = "lm", formula = df$target ~ df$source/(1-df$source^2),col = 'red')+
#  stat_smooth(method = "lm", formula = df$target ~ sinh(df$source),col='green')



#eigenvector evolution
#calculate eigenvector similarity by computing the cosine similarity
sim=cosine(baev_b$vectors[,1],baev_a$vectors[,1])
sim=cosine(baev_b$vectors[,6],baev_a$vectors[,6])

b_svd <- t(svd$u) %*% adj_b 
b_svd<- b_svd %*% svd$v
hist(as.numeric(b_svd))
image(b_svd)
#b_svd <- rescale(b_svd,to = c(0,1))


confusion_matrix <- function(pred,actual,threshold)
{
  if (pred > threshold & actual > 0)
  {
    tp <- tp + 1
  }
  if (pred > threshold & actual == 0)
  {
    fp <- fp + 1
  }
  if (pred < threshold & actual > 0)
  {
    fn <- fn + 1
  }
  if (pred < threshold & actual == 0)
  {
    tn <- tn + 1
  }
  confusion_matrix <- matrix(c(nrow(tp),nrow(fp),nrow(fn),nrow(tn)),nrow = 2,ncol = 2,byrow = TRUE)
  paste(confusion_matrix)
  #  return confusion_matrix
}

apply(mp[1:10,1:10],netest[1:10,1:10],confusion_matrix)

#normalise matrices
normalise <- function(m)
{
  D1 = Matrix(data = 0,nrow=nrow(m),ncol=nrow(m))
  #The function normalises the matrix as described on pages 15-16 
  for (i in 1:nrow(m))
  { 
    #    print(paste('D1 calculation : row', i))
    D1[i,i] = sum(m[i,])
  }
  D2 = Matrix(data = 0,nrow = ncol(m),ncol=ncol(m))  
  for (i in 1:ncol(m))
  {  
    #    print(paste('D2 calculation : column',i))
    D2[i,i] = sum(m[,i])
  }
  diag(D1) <- 1/sqrt(diag(D1))
  diag(D2) <- 1/sqrt(diag(D2))
  M=D1%*%m%*%D2
  return (M)
}

adj_an <- normalise(adj_a)
hist(as.numeric(adj_an))
adj_bn <- normalise(adj_b)
hist(as.numeric(adj_bn))
adj_cn <- normalise(adj_c)
hist(as.numeric(adj_cn))

#spectral evolution test
x<-NULL
y<-NULL
for (i in 1 : 6)
{
  x=append(x,ecount(g_a))
  x=append(x,ecount(g_b))
  x=append(x,ecount(g_c))
  y=append(y,baev_a$values[i])
  y=append(y,baev_b$values[i])
  y=append(y,baev_c$values[i])
}

plot(x,y,xlab = 'Edge count', ylab = 'Eigenvalues')

#compute eigenvector stability
sim_train <- Matrix(nrow = r,ncol = r)
sim_test <- Matrix(nrow = r,ncol = r)
for (i in 1:r)
{
  for (j in 1:r)
  {
    sim_train[i,j] <- abs(cosine(baev_b$vectors[,i],baev_a$vectors[,j]))
    sim_test[i,j] <- abs(cosine(baev_c$vectors[,i],baev_b$vectors[,j]))    
  }
}

image(sim_train,axes=TRUE,main = "Eigenvector stability")
image(sim_test,main = "Eigenvector stability")


#spectral diagonality test
sd <- t(baev_a$vectors) %*% ne %*% baev_a$vectors
image(sd)

#Compare the prediction and the real matrix
ind <- which(mp < 1.304372e-06)
for (i in 1:length(ind))
{
  k <- arrayInd(ind[i],dim(mp))
  rownames(mp)[k[,1]]
  colnames(mp)[k[,2]]
  print(paste('Predicted positive score',mp[ind[i]],' for the the following link'))
  print(mapply(`[[`, dimnames(mp), k))
  print(paste('In reality the following value was observed in the adjacency matrix'))
  print(adj_c[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
  print(paste('New edges matrix had the following value'))
  print(netest[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
}

ind <- which(mp < -0.25)
for (i in 1:length(ind))
{
  k <- arrayInd(ind[i],dim(mp))
  rownames(mp)[k[,1]]
  colnames(mp)[k[,2]]
  print(paste('Predicted negative score',mp[ind[i]],' for the the following link'))
  print(mapply(`[[`, dimnames(mp), k))
  print(paste('In reality the following value was observed in the adjacency matrix'))
  print(adj_c[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
  print(paste('New edges matrix had the following value'))
  print(netest[rownames(mp)[k[,1]],colnames(mp)[k[,2]]])
  print("")
  #  if ((abs(mp[ind[i]])>0) && adj_c[rownames(mp)[k[,1]],colnames(mp)[k[,2]]]>0)
  #    tp <- tp + 1  
}


ppa <- function(m)
  #Probably absolete
  #The function accepts degree matrix and returns rankings based on 
  # preferential attachment model
{
  d <- Matrix(0,nrow = nrow(m),ncol = nrow(m))
  for (i in 1:nrow(d))
  {
    print(paste('Calculating row',i))  
    for (j in 1:ncol(m))
    {  
      if (i!=j) 
      {
        print(paste(i,j))
        d[i,j] = m[i,i]*m[j,j]
        print(d[i])
        print(paste(d[i,j]))
      }
      else 
      {
        d[i,j] = m[i,i]
      }
    }
  }
  return(d)
}
