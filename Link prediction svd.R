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
require(RJDBC)

rm(list=ls(all=TRUE))

jar1 <- "//home//esavin//teradata driver//terajdbc4.jar"
jar2 <- "//home//esavin//teradata driver//tdgssconfig.jar"
jars <- paste(jar1, jar2, sep = ":")

print("opening db connection")
tdconn <- function(id = 'esavin', pwd = 'Lisbon2015!', jarfiles = jars) {
  drv = JDBC("com.teradata.jdbc.TeraDriver", jarfiles)
  conn = dbConnect(drv, "jdbc:teradata://jaguar2.vip.paypal.com/TMODE=TERA", id, pwd)
}

conn = tdconn()

queryS <- paste("select * from pp_oap_sing_es_t.transactions_source")
edgelist_source <- dbGetQuery(conn, queryS)

queryT <- paste("select * from pp_oap_sing_es_t.transactions_target")
edgelist_target <- dbGetQuery(conn, queryT)

queryP <- paste("select * from pp_oap_sing_es_t.transactions_pred")
edgelist_pred <- dbGetQuery(conn, queryP)

print("closing db connection")
dbDisconnect(conn)

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_after.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\ltest_pred.txt',sep = "\t",header = TRUE, row.names = NULL, col.names = c('client','merchant'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_2m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_6m.sorted.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))

#edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_12m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_14m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))
#edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_16m.txt',sep = "|",header = TRUE, row.names = NULL, col.names = c('client','merchant','month'))

#Detect new edges between the target and the source data sets
#new_edges_ts <-edgelist_target[!duplicated(rbind(edgelist_source, edgelist_target))[-seq_len(nrow(edgelist_source))], ]

#Filter out clients active in the source dataset 
#active_clients_source <-edgelist_source[edgelist_source$client %in% new_edges_ts$client,]

#Detect new edges between the prediction and target
#new_edges_pt <-edgelist_pred[!duplicated(rbind(edgelist_target, edgelist_pred))[-seq_len(nrow(edgelist_target))], ]

#Filter out clients active in the source dataset
#active_clients_target <-edgelist_target[edgelist_target$client %in% new_edges_pt$client,]


#sample data
sample_a<-edgelist_source[,c(2,3)] #after data frame
colnames(sample_a) <- c("merchant","client")
sample_b<-edgelist_target[,2:3] #before data frame
colnames(sample_b) <- c("merchant","client")
sample_c<-edgelist_pred[,2:3] #predict data frame
colnames(sample_c) <- c("merchant","client")



#fill the adjacency matrix. 
#clients rows
#merchants columns
g_a <- graph.data.frame(sample_a,directed = FALSE)

adj_a <- get.adjacency(g_a, sparse = TRUE,type=c("upper"))

g_b <- graph.data.frame(sample_b,directed = FALSE)

adj_b <- get.adjacency(g_b, sparse = TRUE,type=c("upper"))

col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_b <- adj_b[row.order,col.order]

g_c <- graph.data.frame(sample_c,directed = FALSE)

adj_c <- get.adjacency(g_c, sparse = TRUE,type=c("upper"))

col.order <- dimnames(adj_a)[[1]]
row.order <- dimnames(adj_a)[[2]]
adj_c <- adj_c[row.order,col.order]



#set the number of eigenvalues to be caluclated
s=100

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

df = as.data.frame(svd_a$d) #the independent variable
df$flamda = as.numeric(diag(b_svd)) #the dependent function
colnames(df) = c("source","target")
plot(df$source,df$target)


#fit <-lm(formula = df$target ~ I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9))

a_poly_odd <-nls(formula = target ~ I(alpha*source) + I(beta*source^3) + I(gamma*source^5),start = c(alpha=0.01,beta=0.01,gamma=0.01), data = df)
summary(a_poly_odd)
sum(resid(a_poly_odd)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_poly_odd,new)
points(x,y,col = "green")

#Still not sure which regression method to use

a_sinh <- nls(formula = target ~ I(alpha*source)+I(((alpha*source)^3)/6)+I(((alpha*source)^5)/120),start=c(alpha=0.01),data=df)
summary(a_sinh)
sum(resid(a_sinh)^2)
confint(a_sinh)
new = data.frame(source = df$source)
x=df$source
y=predict(a_sinh,new)
points(x,y,col = "red")


a_neu_o <- nls(target ~ I(alpha*source) + I((alpha^source)^3) + I((alpha*source)^5),start=c(alpha=0.01),data = df)
summary(a_neu_o)
sum(resid(a_neu_o)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_neu_o,new)
points(x,y,col = "blue")


#Now construct a prediction based on the best fit
flamda_a_poly_odd = matrix(0,  ncol = s,nrow = s)
new.dfb <- data.frame(source = as.numeric(svd_b$d))
diag(flamda_a_poly_odd) = predict(a_poly_odd,new.dfb)
mp_a_poly_odd <- svd_b$u %*% flamda_a_poly_odd %*% t(svd_b$v)
dimnames(mp_a_poly_odd) <- dimnames(adj_a)
mp_a_poly_odd[1:10,1:10]

flamda_a_sinh = matrix(0,  ncol = s,nrow = s)
new.dfb <- data.frame(source = as.numeric(svd_b$d))
diag(flamda_a_sinh) = predict(a_sinh,new.dfb)
mp_a_sinh <- svd_b$u %*% flamda_a_sinh %*% t(svd_b$v)
dimnames(mp_a_sinh) <- dimnames(adj_a)
mp_a_sinh[1:10,1:10]

flamda_a_neu_o = matrix(0,  ncol = s,nrow = s)
new.dfb <- data.frame(source = as.numeric(svd_b$d))
diag(flamda_a_neu_o) = predict(a_neu_o,new.dfb)
mp_a_neu_o <- svd_b$u %*% flamda_a_neu_o %*% t(svd_b$v)
dimnames(mp_a_neu_o) <- dimnames(adj_a)
mp_a_neu_o[1:10,1:10]

netest <- adj_c - adj_b
#hist(as.numeric(netest))
#image(netest) #new edges

netest[netest >= 1] <- 1
netest[netest < 0] <- 0
hist(as.numeric(netest))


pred_a_poly_odd <- prediction(as.vector(mp_a_poly_odd),as.vector(netest))
perf_a_poly_odd <- performance(pred_a_poly_odd, "tpr", "fpr")
plot(perf_a_poly_odd)
precision_recall_a_poly_odd <- performance(pred_a_poly_odd, "prec", "rec")
plot(precision_recall_a_poly_odd)
sensitivity_specificity_a_poly_odd <- performance(pred_a_poly_odd,"sens","spec")
plot(sensitivity_specificity_a_poly_odd)
lift_a_poly_odd <- performance(pred,"lift","rpp")
plot(lift_a_poly_odd)

acc_pred_a_poly_odd <- performance(pred_a_poly_odd,'acc')
f_a_poly_odd <- performance(pred_a_poly_odd,'f')
plot(f_a_poly_odd)
auc_a_poly_odd <- performance(pred_a_poly_odd,"auc")


pred_a_sinh <- prediction(as.vector(mp_a_sinh),as.vector(netest))
perf_a_sinh <- performance(pred_a_sinh, "tpr", "fpr")
plot(perf_a_sinh)
precision_recall_a_sinh <- performance(pred_a_sinh, "prec", "rec")
plot(precision_recall_a_sinh)
sensitivity_specificity_a_sinh <- performance(pred_a_sinh,"sens","spec")
plot(sensitivity_specificity_a_sinh)
lift_a_sinh <- performance(pred,"lift","rpp")
plot(lift_a_sinh)

acc_pred_a_sinh <- performance(pred_a_sinh,'acc')
err_pred_a_sinh <- performance(pred_a_sinh,'err')
f_a_sinh <- performance(pred_a_sinh,'f')
plot(f_a_sinh)
auc_a_sinh <- performance(pred_a_sinh,"auc")


pred_a_neu_o <- prediction(as.vector(mp_a_neu_o),as.vector(netest))
perf_a_neu_o <- performance(pred_a_neu_o, "tpr", "fpr")
plot(perf_a_neu_o)
precision_recall_a_neu_o <- performance(pred_a_neu_o, "prec", "rec")
plot(precision_recall_a_neu_o)
sensitivity_specificity_a_neu_o <- performance(pred_a_neu_o,"sens","spec")
plot(sensitivity_specificity_a_neu_o)
lift_a_neu_o <- performance(pred,"lift","rpp")
plot(lift_a_neu_o)

acc_pred_a_neu_o <- performance(pred_a_neu_o,'acc')
f_a_neu_o <- performance(pred_a_neu_o,'f')
plot(f_a_neu_o)
auc_a_neu_o <- performance(pred_a_neu_o,"auc")

cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
                      tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
head(cutoffs)
head(subset(cutoffs, fpr < 0.50))

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

#calculate preferential attachment scores
adj_b_d <- degree(adj_b)

adj_b_ppa<-as.matrix(as.vector(diag(adj_b_d))%*%t(as.vector(diag(adj_b_d))))

h<- hist(adj_b_ppa)
h$counts

pred_ppa <- prediction(as.vector(adj_b_ppa),as.vector(netest))
perf_ppa <- performance(pred_ppa,"tpr","fpr")
plot(perf_ppa)
sensitivity_specificity_ppa <- performance(pred,"sens","spec")
plot(sensitivity_specificity_ppa)
precision_recall_ppa <- performance(pred_ppa, "prec", "rec")
plot(precision_recall_ppa)
auc_ppa <- performance(pred_ppa,"auc")
auc_ppa


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
    print(paste('D2 calculation : column',i))
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


############################################################################



#find concetration of values and set them to be 0
#mp<-abs(mp)
h<-hist(as.numeric(mp))
h$counts
#z<-match(max(h$counts),h$counts) #zero 
#h$breaks[z]
l<-length(h$breaks)
#mp[mp<h$breaks[z-l]] <- 1
mp[mp<=h$breaks[l-round(l/1.2)]] <- 0
mp[mp>h$breaks[l-round(l/1.2)]] <- 1
mp[mp<=0.0006290906] <- 0
mp[mp>0.0006290906] <- 1
hist(as.numeric(mp))
#image(mp)
#image(netest)


ind <- which(mp == 1)
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


ind <- which(mp > 0.5)
#substitute any multiple edges with just 1 edge

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

tpr = nrow(tp)/(nrow(tp)+nrow(fp))
tnr = nrow(tn)/(nrow(tn)+nrow(fn))
fpr = nrow(fp)/(nrow(fp)+nrow(tn))
fnr = nrow(fn)/(nrow(fn)+nrow(tn))

precision = nrow(tp)/(nrow(tp)+nrow(fp)) #positive predictive value
npv = nrow(tn)/(nrow(fn)+nrow(tn)) #negative predictive value
recall = nrow(tp)/(nrow(fn)+nrow(tp)) #sensitivity
specificity = nrow(tn)/(nrow(fp)+nrow(tn))

precision_matrix <- matrix(c(precision,npv,recall,specificity),nrow = 2,ncol = 2,byrow = TRUE)
precision_matrix

f_measure_pos = (2*precision*tpr)/(precision+tpr)
f_measure_neg = (2*specificity*tnr)/(specificity+tnr)


g_mean = sqrt(tpr*tnr)




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


#detect abnormally active merchants and filter them out
a <- as.data.frame(table(sample_a$merchant),row.names = NULL)
a<- a[order(a$Freq,decreasing = TRUE),]
hist(as.numeric(a$Freq))

for (i in 1:nrow(a))
  if (a$Freq[i] > (mean(a$Freq) + 3*sd(a$Freq)) || (a$Freq[i] < 5))
  {
    #print(a$Var1[i])
    sample_a <- sample_a[sample_a$merchant != a$Var1[i],]
    sample_b <- sample_b[sample_b$merchant != a$Var1[i],]
    sample_c <- sample_c[sample_b$merchant != a$Var1[i],]  
  }


outersect <- function (x,y)
{
  sort(c(setdiff(x,y),setdiff(y,x)))
}


diff_as1<-outersect(edgelist_source$counterparty_id,edgelist_target$counterparty_id)
diff_as2<-outersect(edgelist_source$counterparty_id,edgelist_pred$counterparty_id)


diff_as1<-outersect(edgelist_source$cust_id,edgelist_target$cust_id)
diff_as2<-outersect(edgelist_source$cust_id,edgelist_pred$cust_id)


setdiff(edgelist_source$counterparty_id,edgelist_target$counterparty_id)
setdiff(edgelist_target$counterparty_id,edgelist_source$counterparty_id)
#Equate the number of nodes
#Run a number of times until diff1 and diff2 are empy
diff1_as1 <- outersect(sample_b$client,sample_a$client)
diff2_as2 <- outersect(sample_b$merchant,sample_a$merchant)

#Remove rows where elements appear in just 1 data set and not both, i.e. make sure both datasets
#have the same number of nodes
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


a_poly_odd <-lm(formula = target ~ I(source) + I(source^3) + I(source^5), data = df)
summary(a_poly_odd)
sum(resid(a_poly_odd)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(a_poly_odd,new)
points(x,y,col = "green")

alpha = 0.1
a_sinh <- lm(formula = target ~ I(alpha*source)+I(((alpha*source)^3)/6)+I(((alpha*source)^5)/120),data=df)
summary(a_sinh)
sum(resid(a_sinh)^2)
confint(a_sinh)
new = data.frame(source = df$source)
x=df$source
y=predict(a_sinh,new)
points(x,y,col = "red")

a_neu_o <- nls(formula = target ~ I(alpha*source) + I(alpha^3*source^3) + I(alpha^5*source^5),start=c(alpha=0.01),data = df)
summary(a_neu_o)
sum(resid(a_neu_o)^2)
confint(a_neu_o)
new = data.frame(xdata = seq(min(df$source),max(df$source),len=length(df$source)))
x=df$source
y=predict(a_neu_o,new)
points(x,y,col = "blue")


#adj_a <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
#rownames(adj_a) <- c(as.character(unique(sample_a$client)))
#colnames(adj_a) <- c(as.character(unique(sample_a$merchant)))

#for (i in 1:nrow(sample_a))
#  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] <- 1
#  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] <- 
#  adj_a[as.character(sample_a$client[i]),as.character(sample_a$merchant[i])] +1

#image(adj_a)
#hist(as.numeric(adj_a))

#adj_b <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
#rownames(adj_b) <- c(as.character(unique(sample_a$client)))
#colnames(adj_b) <- c(as.character(unique(sample_a$merchant)))

#for (i in 1:nrow(sample_b))
#{
#unweighted graph
#  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 1
#weighted graph 
#  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] <- 
#  adj_b[as.character(sample_b$client[i]),as.character(sample_b$merchant[i])] +1
#}

#adj_c <- get.adjacency(g_c, sparse = TRUE,type=c("upper"))

#adj_c <- Matrix(0,sparse = TRUE,nrow = length(unique(sample_a$client)), ncol = length(unique(sample_a$merchant)))
#rownames(adj_c) <- c(as.character(unique(sample_a$client)))
#colnames(adj_c) <- c(as.character(unique(sample_a$merchant)))

#for (i in 1:nrow(sample_c))
#unweighted graph
#  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] <- 1
#weighted graph 
#  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] <- 
#  adj_c[as.character(sample_c$client[i]),as.character(sample_c$merchant[i])] +1
