require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)
require(SnowballC)
require(lsa)
require(scales)
library(reshape2)
options(java.parameters = "-Xmx4g")

rm(list=ls(all=TRUE))

edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_1.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_3.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_4m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_after.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_pred.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_12m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_14m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_16m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))


#edgelist_source<-edgelist_source[order(edgelist_source$client),]
#edgelist_target<-edgelist_target[order(edgelist_target$client),]
#edgelist_pred<-edgelist_pred[order(edgelist_pred$client),]


#Detect new edges between the target and the source data sets
new_edges_ts <-edgelist_target[!duplicated(rbind(edgelist_source, edgelist_target))[-seq_len(nrow(edgelist_source))], ]
new_edges_train <- as.data.frame(new_edges_ts$client) 
new_edges_train$merchant <- new_edges_ts$merchant
colnames(new_edges_train) <- c("client","merchant")
new_edges_ts <- NULL

#Filter out clients active in the source dataset 
#active_clients_source <-edgelist_source[edgelist_source$client %in% new_edges_ts$client,]

#Detect new edges between the prediction and target
new_edges_pt <-edgelist_pred[!duplicated(rbind(edgelist_target, edgelist_pred))[-seq_len(nrow(edgelist_target))], ]
new_edges_test <- as.data.frame(new_edges_pt$client) 
new_edges_test$merchant <- new_edges_pt$merchant
new_edges_pt <- NULL
colnames(new_edges_test) <- c("client","merchant")

#Filter out clients active in the source dataset
#active_clients_target <-edgelist_target[edgelist_target$client %in% new_edges_pt$client,]


#sample data
sample_a<-edgelist_source #after data frame
sample_b<-edgelist_target #before data frame
sample_c<-edgelist_pred #predict data frame


#set the number of eigenvalues to be caluclated
r=100

g_a <- graph.data.frame(sample_a,directed = FALSE)
V(g_a)$type <- V(g_a)$name %in% sample_a[,1]
is.bipartite(g_a) #verify graph is bipartite
get.edgelist(g_a, names=TRUE)[1:10,]# to verify this was done properly
adj_a <- get.adjacency(g_a, sparse = TRUE,type=c("both"))
image(adj_a)


#isSymmetric(as.matrix(adj))
#calculate the eigenvalues and eigenvectors using ARPACK algorithm

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a %*% x) }
baev_a <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))


g_b <- graph.data.frame(sample_b,directed = FALSE)
V(g_b)$type <- V(g_b)$name %in% sample_b[,1]
#get.edgelist(g_b) #to verify this was done properly
adj_b <- get.adjacency(g_b, sparse = TRUE)
#reorder matrix adj_b to match order of adj_a
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_b <- adj_b[row.order,col.order]
image(adj_b)

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b %*% x) }
baev_b <- arpack(f2, sym=TRUE, options=list(n=vcount(g_b), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_b)*12))


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
V(g_c)$type <- V(g_c)$name %in% sample_c[,1]
#get.edgelist(g_c) #to verify this was done properly
adj_c <- get.adjacency(g_c, sparse = TRUE)
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_c <- adj_c[row.order,col.order]
f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_c %*% x) }
baev_c <- arpack(f2, sym=TRUE, options=list(n=vcount(g_c), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_c)*15))
mc = Matrix(0,  ncol = r,nrow = r)
diag(mc) = baev_c$values


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


#compute multiplication of eigenvalues of A and matrix B
b_eigen <- t(baev_a$vectors) %*% adj_b %*% baev_a$vectors
#hist(as.numeric(b_eigen))
image(b_eigen)


#df = as.data.frame(nlamda)
df = as.data.frame(as.numeric(diag(b_eigen))) #the independent variable
df$flamda = nlamda #the dependent function
colnames(df) = c("source","target")
plot(df$source,df$target)


#fit1 <-lm(formula = df$target ~ I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9))
fit1 <-lm(formula = target ~ I(source) + I(source^3) + I(source^5) + I(source^7) + I(source^9), data = df)
summary(fit1)
sum(resid(fit1)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(fit1,new)
points(x,y,col = "green")


fit2 <-lm(formula = target ~ sinh(I(source)),data = df)
summary(fit2)
sum(resid(fit2)^2)
confint(fit2)
new = data.frame(source = df$source)
x=df$source
y=predict(fit2,new)
points(x,y,col = "red")

fit3 <- lm(formula = target ~ I(source)/(1-I(source^2)), data = df)
summary(fit3)
sum(resid(fit3)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(fit3,new)
lines(x,y,col = "blue")


netest <- adj_c - adj_b
image(netest) #new edges


#Now construct a prediction based on the best fit
flamda = Matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
newdata = predict(fit1,new.dfb)
#newdata = coef(fit)[1] + coef(fit)[2]*baev_b$values + (baev_b$values^3)*coef(fit)[3] + coef(fit)[4]*(baev_b$values^5) + coef(fit)[5]*(baev_b$values^7) + coef(fit)[6]*(baev_b$values^9)
hist(newdata)
hist(as.numeric(mc))
diag(flamda) = newdata
image(mc)
image(flamda)
mp <- baev_b$vectors %*% flamda %*% t(baev_b$vectors)
mp<-abs(mp)

dimnames(mp) <- dimnames(adj_a)
h<-hist(as.numeric(mp))
#find concetration of values and set them to be 0
z<-match(max(h$counts),h$counts) #zero 
l<-length(h$counts)
#h$breaks[z]
#mp[mp<=h$breaks[1]] <- 0
mp[mp<=h$breaks[l-4]] <- 0
mp[mp>h$breaks[l-4]] <- 1
hist(as.numeric(mp))
image(mp)
netest <- adj_c - adj_b
image(netest)

m<-netest+2*mp
hist(as.numeric(m))
#m<-2*as.matrix(mp)+as.matrix(netest)
#image(m) #new edges
tp <- which(m == 3, arr.in = TRUE, useNames = TRUE)
fp <- which(m == 2, arr.in = TRUE, useNames = TRUE)
fn <- which(m == 1, arr.in = TRUE, useNames = TRUE)
tn <- which(m == 0, arr.in = TRUE, useNames = TRUE)

confusion_matrix <- matrix(c(nrow(tp),nrow(fp),nrow(fn),nrow(tn)),nrow = 2,ncol = 2,byrow = TRUE)
confusion_matrix


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
preddf$merchant <- pred$Var1
colnames(preddf) <- c("client","merchant")

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
preddf$merchant <- pred$Var1
colnames(preddf) <- c("client","merchant")

match <-new_edges_test[!duplicated(rbind(new_edges_test, preddf))[-seq_len(nrow(new_edges_test))], ]









###############################################################################################


#function for selecting a subsample of transacting customers and merchants from 3 different datasets
sample_a <- edgelist_source
sample_a <- NULL
sample_b <- edgelist_source
sample_b <- NULL
sample_c <- edgelist_source
sample_c <- NULL

x <- as.integer(runif(10,1,nrow(edgelist_source))) 
#x <- active_clients_source[1:5,]
#d<-edgelist_source[1:5,]$client
d<-edgelist_source[x,]$client
for (i in 1 : length(d))
{
  #i=1
  sc_source<-edgelist_source[edgelist_source$client==d[i],] #set of clients source
  sc_target<-edgelist_target[edgelist_target$client==d[i],] #set of clients target
  sc_pred<-edgelist_pred[edgelist_pred$client==d[i],] #set of clients prediction
  
  cst<-setdiff(sc_target$merchant,sc_source$merchant) #check that there were no new merchants between source and target 
  ctp<-setdiff(sc_pred$merchant,sc_source$merchant) #check that there were no new merchants added between target and prediction 
  
  for (j in 1 : length(cst))
  {
    sc_target <- sc_target[sc_target$merchant != cst[j],]
  }
  
  for (l in 1 : length(ctp))
  {
    sc_pred <- sc_pred[sc_pred$merchant != ctp[l],]
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


#verify the merchant sets are identical
setdiff(sample_b$merchant,sample_a$merchant)
setdiff(sample_c$merchant,sample_b$merchant)

setdiff(sample_a$merchant,sample_b$merchant)
setdiff(sample_b$merchant,sample_c$merchant)

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

#singular value decomposition
require(irlba)
s=10
svd<-irlba(adj_a, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
           sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
           tol = 1e-06, V = NULL, matmul = NULL)
svd_m <- Matrix(0, ncol = s,nrow = s)
diag(svd_m) <- svd$d
t<-svd$u %*% svd_m 
t<-t %*% t(svd$v)
image(adj_a)
t <- rescale(t,to = c(0,1))
hist(as.numeric(t))
t[t<0.4] <- 0
t[t>0.4] <- 1
t <- Matrix(t,sparse = TRUE)
image(t)


#eigenvector evolution
#calculate eigenvector similarity by computing the cosine similarity
sim=cosine(baev_b$vectors[,1],baev_a$vectors[,1])
sim=cosine(baev_b$vectors[,6],baev_a$vectors[,6])

b_svd <- t(svd$u) %*% adj_b 
b_svd<- b_svd %*% svd$v
hist(as.numeric(b_svd))
image(b_svd)
#b_svd <- rescale(b_svd,to = c(0,1))

