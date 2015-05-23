require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)
require(SnowballC)
require(lsa)
require(scales)
library(reshape2)
require(irlba)
options(java.parameters = "-Xmx4g")

rm(list=ls(all=TRUE))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_after.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_pred.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))

edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_12m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_14m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_16m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))

#Detect new edges between the target and the source data sets
new_edges_ts <-edgelist_target[!duplicated(rbind(edgelist_source, edgelist_target))[-seq_len(nrow(edgelist_source))], ]

#Filter out clients active in the source dataset 
#active_clients_source <-edgelist_source[edgelist_source$client %in% new_edges_ts$client,]

#Detect new edges between the prediction and target
new_edges_pt <-edgelist_pred[!duplicated(rbind(edgelist_target, edgelist_pred))[-seq_len(nrow(edgelist_target))], ]

#Filter out clients active in the source dataset
#active_clients_target <-edgelist_target[edgelist_target$client %in% new_edges_pt$client,]


#sample data
sample_a<-edgelist_source[,1:2] #after data frame
sample_b<-edgelist_target[,1:2] #before data frame
sample_c<-edgelist_pred[,1:2] #predict data frame

outersect <- function (x,y)
{
  sort(c(setdiff(x,y),setdiff(y,x)))
}


#detect abnormally active clients and filter them out
a <- as.data.frame(table(sample_a$client),row.names = NULL)
a<- a[order(a$Freq,decreasing = TRUE),]
hist(as.numeric(a$Freq))

for (i in 1:nrow(a))
{
  if (a$Freq[i] > 60 || a$Freq[i] < 4)
  {
    #print(a$Var1[i])
    sample_a <- sample_a[sample_a$client != a$Var1[i],]
    sample_b <- sample_b[sample_b$client != a$Var1[i],]
    sample_c <- sample_c[sample_c$client != a$Var1[i],]  
  }
}  

outersect(sample_b$client,sample_a$client)
outersect(sample_b$merchant,sample_a$merchant)


#detect abnormally active merchants and filter them out
a <- as.data.frame(table(sample_a$merchant),row.names = NULL)
a<- a[order(a$Freq,decreasing = TRUE),]
hist(as.numeric(a$Freq))

for (i in 1:nrow(a))
  if (a$Freq[i] > (mean(a$Freq) + sd(a$Freq)) || (a$Freq[i] < (mean(a$Freq) - sd(a$Freq))))
  {
    #print(a$Var1[i])
    sample_a <- sample_a[sample_a$merchant != a$Var1[i],]
    sample_b <- sample_b[sample_b$merchant != a$Var1[i],]
    sample_c <- sample_c[sample_b$merchant != a$Var1[i],]  
  }

#Equate the number of nodes
#Run a number of times until diff1 and diff2 are empy
diff1 <- outersect(sample_b$client,sample_a$client)
diff2 <- outersect(sample_b$merchant,sample_a$merchant)

#Remove rows where elements appear in just 1 data set and not both, i.e. make sure both datasets
#have the same number of nodes
for (i in 1:length(diff1))
{
  sample_a <- sample_a[sample_a$client != diff1[i],]
  sample_b <- sample_b[sample_b$client != diff1[i],]
}

for (i in 1:length(diff2))
{
  sample_a <- sample_a[sample_a$merchant != diff2[i],]
  sample_b <- sample_b[sample_b$merchant != diff2[i],]
}

diff1 <- outersect(sample_c$client,sample_a$client)
diff2 <- outersect(sample_c$merchant,sample_a$merchant)

#Remove rows where elements appear in just 1 data set and not both, i.e. make sure both datasets
#have the same number of nodes
for (i in 1:length(diff1))
{
  sample_a <- sample_a[sample_a$client != diff1[i],]
  sample_c <- sample_c[sample_c$client != diff1[i],]
}

for (i in 1:length(diff2))
{
  sample_a <- sample_a[sample_a$merchant != diff2[i],]
  sample_c <- sample_c[sample_c$merchant != diff2[i],]
}



#set the number of eigenvalues to be caluclated
s=80

#fill the adjacency matrix. 
#clients rows
#merchants columns
adj_a <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
rownames(adj_a) <- c(as.character(unique(sample_a$client)))
colnames(adj_a) <- c(as.character(unique(sample_a$merchant)))

for (i in 1:nrow(sample_a))
#  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] <- 1
  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] <- 
  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] +1

#image(adj_a)
hist(as.numeric(adj_a))

adj_b <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
rownames(adj_b) <- c(as.character(unique(sample_a$client)))
colnames(adj_b) <- c(as.character(unique(sample_a$merchant)))

for (i in 1:nrow(sample_b))
{
  #unweighted graph
#  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 1
#weighted graph 
  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 
  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] +1
}

#image(adj_b)
hist(as.numeric(adj_b))

adj_c <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
rownames(adj_c) <- c(as.character(unique(sample_a$client)))
colnames(adj_c) <- c(as.character(unique(sample_a$merchant)))

for (i in 1:nrow(sample_c))
#unweighted graph
#  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] <- 1
#weighted graph 
  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] <- 
  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] +1

#image(adj_c)
hist(as.numeric(adj_c))


#matrix normalisation
normalise <- function(m)
{
  D1 = Matrix(data = 0,nrow=nrow(m),ncol=nrow(m))
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
adj_bn <- normalise(adj_b)
hist(as.numeric(adj_bn))
adj_cn <- normalise(adj_c)
hist(as.numeric(adj_cn))


svd_a<-irlba(adj_a, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
#svd_am <- Matrix(0, ncol = s,nrow = s)
#diag(svd_am) <- svd_a$d
#t<-svd_a$u %*% svd_am %*% t(svd_a$v)
#image(adj_a)
##image(t)

svd_b<-irlba(adj_b, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
#svd_bm <- Matrix(0, ncol = s,nrow = s)
#diag(svd_bm) <- svd_b$d
#t<-svd_b$u %*% svd_bm %*% t(svd_b$v)
#image(adj_b)
##image(t)

#spectral evolution test
plot(x=c(nrow(sample_a),nrow(sample_b)),y=c(svd_a$d[1],svd_b$d[1]),xlab = 'Edge count ', ylab = 'Singular values')

ne <- adj_b - adj_a
#image(adj_a)
#image(adj_b)
#image(ne)

svd_ne<-irlba(ne, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
#svd_nem <- Matrix(0, ncol = s,nrow = s)
#diag(svd_nem) <- svd_ne$d
#t<-svd_ne$u %*% svd_nem %*% t(svd_ne$v)
##image(ne)
##image(t)


#spectral diagonality test
sd <- t(svd_a$u) %*% ne %*% svd_a$v
dim(sd)
image(sd)

#normalise the spectra
nlamda <- svd_a$d * abs(svd_a$d[1])/abs(svd_b$d[1])

#without normalisation
#nlamda <- svd_a$d

#compute multiplication of eigenvalues of A and matrix of new edges of B
#ne <- adj_b
b_svd <- t(svd_a$u) %*% ne %*% svd_a$v
#hist(as.numeric(b_svd))
image(b_svd)

df = as.data.frame(as.numeric(diag(b_svd))) #the independent variable
df$flamda = nlamda #the dependent function
colnames(df) = c("source","target")
plot(df$source,df$target)


#fit <-lm(formula = df$target ~ I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9))
fit1 <-lm(formula = target ~ I(source) + I(source^3) + I(source^5) + I(source^7) + I(source^9), data = df)
summary(fit1)
sum(resid(fit1)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(fit1,new)
points(x,y,col = "green")


fit2 <- lm(formula = df$target ~ sinh(df$source))
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
points(x,y,col = "blue")


netest <- adj_c - adj_b
#image(netest) #new edges

#Now construct a prediction based on the best fit
flamda = Matrix(0,  ncol = s,nrow = s)
new.dfb <- data.frame(source = as.numeric(svd_b$d))
newdata = predict(fit1,new.dfb)
#newdata = coef(fit)[1] + coef(fit)[2]*baev_b$values + (baev_b$values^3)*coef(fit)[3] + coef(fit)[4]*(baev_b$values^5) + coef(fit)[5]*(baev_b$values^7) + coef(fit)[6]*(baev_b$values^9)
#hist(newdata)
diag(flamda) = newdata
##image(svd_cm)
##image(flamda)
mp <- svd_b$u %*% flamda %*% t(svd_b$v)
dimnames(mp) <- dimnames(adj_a)
#find concetration of values and set them to be 0
mp<-abs(mp)
h<-hist(as.numeric(mp))
#z<-match(max(h$counts),h$counts) #zero 
#h$breaks[z]
l<-length(h$breaks)
#mp[mp<h$breaks[z-l]] <- 1
mp[mp<=h$breaks[l-round(l/1.2)]] <- 0
mp[mp>h$breaks[l-round(l/1.2)]] <- 1
hist(as.numeric(mp))
#image(mp)
#image(netest)

ind <- which(mp == 1)
for (i in 1:length(ind))
{
  k <- arrayInd(ind[i],dim(mp))
  print(mapply(`[[`, dimnames(mp), k))
}

#substitute any multiple edges with just 1 edge
netest[netest > 1] <- 1
mx<-max(as.numeric(netest))

m<-netest+1.5*mp
hist(as.numeric(m))
#m<-2*as.matrix(mp)+as.matrix(netest)
##image(m) #new edges

tp <- which(m == mx+1.5, arr.in = TRUE, useNames = TRUE)
fp <- which(m == 1.5, arr.in = TRUE, useNames = TRUE)
fn <- which(m == mx, arr.in = TRUE, useNames = TRUE)
tn <- which(m == 0, arr.in = TRUE, useNames = TRUE)

confusion_matrix <- matrix(c(nrow(tp),nrow(fp),nrow(fn),nrow(tn)),nrow = 2,ncol = 2,byrow = TRUE)
confusion_matrix

precision = nrow(tp)/(nrow(tp)+nrow(fp)) #positive predictive value
npv = nrow(tn)/(nrow(fn)+nrow(tn)) #negative predictive value
recall = nrow(tp)/(nrow(fn)+nrow(tp)) #sensitivity
specificity = nrow(tn)/(nrow(fp)+nrow(tn))

precision_matrix <- matrix(c(precision,npv,recall,specificity),nrow = 2,ncol = 2,byrow = TRUE)
precision_matrix

pred <- melt(as.matrix(mp))
p1 <- pred[pred$value == 1,]



