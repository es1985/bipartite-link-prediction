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
edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_4m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_after.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_pred.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

edgelist_source<-edgelist_source[order(edgelist_source$client),]
edgelist_target<-edgelist_target[order(edgelist_target$client),]
edgelist_pred<-edgelist_pred[order(edgelist_pred$client),]


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
s=300

g_a <- graph.data.frame(sample_a,directed = FALSE)
V(g_a)$type <- V(g_a)$name %in% sample_a[,1]
#is.bipartite(g_a) #verify graph is bipartite
#get.edgelist(g, names=TRUE)# to verify this was done properly
adj_a <- get.adjacency(g_a, sparse = TRUE,type=c("both"))
#adj_a[is.na(adj_a)] <- 0
#image(adj_a,zlim = c(0,1))
#nrow(adj_a)
#ncol(adj_a)

svd_a<-irlba(adj_a, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_am <- Matrix(0, ncol = s,nrow = s)
diag(svd_am) <- svd_a$d
#t<-svd_a$u %*% svd_am %*% t(svd_a$v)
#hist(as.numeric(t))
#t <- Matrix(t,sparse = TRUE)
#image(adj_a)
#image(t, zlim = c(0,1))

#isSymmetric(as.matrix(adj))
#calculate the eigenvalues and eigenvectors using ARPACK algorithm

g_b <- graph.data.frame(sample_b,directed = FALSE)
V(g_b)$type <- V(g_b)$name %in% sample_b[,1]
#get.edgelist(g_b) #to verify this was done properly
adj_b <- get.adjacency(g_b, sparse = TRUE,type=c("both"))
#reorder matrix adj_b to match order of adj_a
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_b <- adj_b[row.order,col.order]
#adj_b[is.na(adj_b)] <- 0
#image(adj_b)

#compute matrix of new edges - i.e. B
#ne <- as.matrix(adj_b) - as.matrix(adj_a)
ne <- adj_b - adj_a
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
ne <- ne[row.order,col.order]
#ne <- Matrix(ne,sparse = TRUE)
#image(ne,xaxs = "i",yaxs = "i")
#image(adj_a)
#image(adj_b)
#image(ne)
#adj_b <- ne
#adj_b <- ne
#isSymmetric(as.matrix(adj))
#calculate the eigenvalues and eigenvectors using ARPACK algorithm
svd_b<-irlba(adj_b, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_bm <- Matrix(0, ncol = s,nrow = s)
diag(svd_bm) <- svd_b$d
#t<-svd_b$u %*% svd_bm %*% t(svd_b$v)
#hist(as.numeric(t))
#image(adj_b)
#image(t)

#compute eigenvalues for the new edges matrix of B.
svd_ne<-irlba(ne, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
              sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
              tol = 1e-06, V = NULL, matmul = NULL)
svd_nem <- Matrix(0, ncol = s,nrow = s)
diag(svd_nem) <- svd_ne$d
#t<-svd_ne$u %*% svd_nem %*% t(svd_ne$v)
#hist(as.numeric(t))
#image(adj_b)
#image(ne)
#image(t)

#spectral evolution test
plot(x=c(ecount(g_a),ecount(g_b)),y=c(svd_a$d[1],svd_b$d[1]),xlab = 'Eigenvalues', ylab = 'Edge count')

#spectral diagonality test
sd <- t(svd_a$v) %*% adj_b %*% svd_a$u
dim(sd)
#image(sd)

#Check the development of the prediction eigenvalues in comparison to target
g_c <- graph.data.frame(sample_c,directed = FALSE)
V(g_c)$type <- V(g_c)$name %in% sample_c[,1]
#get.edgelist(g_c) #to verify this was done properly
adj_c <- get.adjacency(g_c, sparse = TRUE,type = "both")
col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_c <- adj_c[row.order,col.order]
svd_c<-irlba(adj_c, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_cm <- Matrix(0, ncol = s,nrow = s)
diag(svd_cm) <- svd_c$d
#t<-svd_c$u %*% svd_cm %*% t(svd_c$v)
#hist(as.numeric(t))
#image(adj_c)
#image(t)

#normalise the spectra
nlamda <- svd_a$d * abs(svd_a$d[1])/abs(svd_b$d[1])

#without normalisation
#nlamda <- svd_a$d

#compute multiplication of eigenvalues of A and matrix of new edges of B
#ne <- adj_b
b_svd <- t(svd_a$v) %*% ne %*% svd_a$u
#hist(as.numeric(b_svd))
#image(b_svd)

df = as.data.frame(as.numeric(diag(b_svd))) #the independent variable
df$flamda = nlamda #the dependent function
colnames(df) = c("source","target")
plot(df$source,df$target)


#fit <-lm(formula = df$target ~ I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9))
fit1 <-lm(formula = target ~ I(source) + I(source^3) + I(source^5) + I(source^7) + I(source^9), data = df)
summary(fit1)
sum(resid(fit)^2)
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
lines(x,y,col = "blue")


#Sanity check - reproduce the target matrix
#flamda = Matrix(0,  ncol = ncol(b_svd),nrow = ncol(b_svd),sparse=TRUE)
#diag(flamda) = coef(fit)[1] + coef(fit)[2]*baev_a$values + (baev_a$values^3)*coef(fit)[3] + coef(fit)[4]*(baev_a$values^5) + coef(fit)[5]*(baev_a$values^7) + coef(fit)[6]*(baev_a$values^9)
#image(mb)
#image(flamda)
#mp <- Matrix(baev_a$vectors[,1:r] %*% flamda %*% t(baev_a$vectors[,1:r]),sparse=TRUE)
#hist(as.numeric(mp))
#substitute values in a matrix
#image(adj_a)
#image(adj_b)
#image(ne)
#image(mp)

netest <- adj_c - adj_b
image(netest) #new edges

#Now construct a prediction based on the best fit
flamda = Matrix(0,  ncol = s,nrow = s)
new.dfb <- data.frame(source = as.numeric(svd_b$d))
newdata = predict(fit1,new.dfb)
#newdata = coef(fit)[1] + coef(fit)[2]*baev_b$values + (baev_b$values^3)*coef(fit)[3] + coef(fit)[4]*(baev_b$values^5) + coef(fit)[5]*(baev_b$values^7) + coef(fit)[6]*(baev_b$values^9)
hist(newdata)
diag(flamda) = newdata
#image(svd_cm)
#image(flamda)
mp <- svd_b$u %*% flamda %*% t(svd_b$v)
#find concetration of values and set them to be 0
mp<-abs(mp)
h<-hist(as.numeric(mp))
z<-match(max(h$counts),h$counts) #zero 
#h$breaks[z]
l<-length(h$breaks)
round(l/2)
#mp[mp<h$breaks[z-l]] <- 1
mp[mp<=h$breaks[l-16]] <- 0
mp[mp>h$breaks[l-16]] <- 1
hist(as.numeric(mp))
image(mp)
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

precision = nrow(tp)/(nrow(tp)+nrow(fp)) #positive predictive value
npv = nrow(tn)/(nrow(fn)+nrow(tn)) #negative predictive value
recall = nrow(tp)/(nrow(fn)+nrow(tp)) #sensitivity
specificity = nrow(tn)/(nrow(fp)+nrow(tn))

precision_matrix <- matrix(c(precision,npv,recall,specificity),nrow = 2,ncol = 2,byrow = TRUE)
precision_matrix

#image(adj_c)
#new edges test

#image(ne)
#image(adj_a)
#image(adj_b)
#image(adj_c)
max(as.numeric(mp))
mean(as.numeric(mp))
min(as.numeric(mp))

mp <- as.matrix(mp,row.names = FALSE)


yo <- which(mp>0.8,arr.ind = TRUE,useNames=FALSE)

pred <- as.matrix(mp[yo[,1],yo[,2]])
hist(as.numeric(pred))

t <- rescale(as.numeric(pred),to = c(0,1))
hist(as.numeric(t))
t[t<0.4] <- 0
t[t>0.4] <- 1

pred <- melt(pred)
pred <- pred[,1:2]
colnames(pred) <- c("client","merchant")

pred <- pred[pred$value == 1,]
pred <- pred[order(pred[1,])]
preddf <- NULL
preddf <-as.data.frame(pred$Var2)
preddf$merchant <- pred$Var1
colnames(preddf) <- c("client","merchant")

dimnames(m==2)

m>0


which(pred>0.8)



image(mp[1:100,1:100]) #predicted new edges
