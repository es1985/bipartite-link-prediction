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

edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_12m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_14m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_16m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

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
sample_a<-edgelist_source[,1:2] #after data frame
sample_b<-edgelist_target[,1:2] #before data frame
sample_c<-edgelist_pred #predict data frame


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

image(adj_a)
hist(as.numeric(adj_a))

adj_b <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
rownames(adj_b) <- c(as.character(unique(sample_a$client)))
colnames(adj_b) <- c(as.character(unique(sample_a$merchant)))

for (i in 1:nrow(sample_b))
#unweighted graph
#  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 1
#weighted graph 
  i=1
  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 
  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] +1

image(adj_b)
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

image(adj_c)
hist(as.numeric(adj_c))

svd_a<-irlba(adj_a, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_am <- Matrix(0, ncol = s,nrow = s)
diag(svd_am) <- svd_a$d
t<-svd_a$u %*% svd_am %*% t(svd_a$v)
image(adj_a)
image(t)

svd_b<-irlba(adj_b, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_bm <- Matrix(0, ncol = s,nrow = s)
diag(svd_bm) <- svd_b$d
t<-svd_b$u %*% svd_bm %*% t(svd_b$v)
image(adj_b)
image(t)

#spectral evolution test
plot(x=c(nrow(sample_a),nrow(sample_b)),y=c(svd_a$d[1],svd_b$d[1]),xlab = 'Edge count ', ylab = 'Singular values')

ne <- adj_b - adj_a
image(adj_a)
image(adj_b)
image(ne)

svd_ne<-irlba(ne, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
             sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
             tol = 1e-06, V = NULL, matmul = NULL)
svd_nem <- Matrix(0, ncol = s,nrow = s)
diag(svd_nem) <- svd_ne$d
t<-svd_ne$u %*% svd_nem %*% t(svd_ne$v)
image(ne)
image(t)


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
lines(x,y,col = "blue")

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
mp[mp<=h$breaks[l-5]] <- 0
mp[mp>h$breaks[l-5]] <- 1
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
