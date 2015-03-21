require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)

rm(list=ls(all=TRUE))
edgelist_source <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_12m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_target <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_14m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))
edgelist_pred <- read.csv('C:\\Users\\esavin\\Documents\\Link prediction\\pp_oap_sing_tv_t.graph_edges_16m.sorted.txt',sep = ",",header = TRUE, row.names = NULL, col.names = c('client','merchant'))


#Detect new edges between the target and the source data sets
new_edges_ts <-edgelist_target[!duplicated(rbind(edgelist_source, edgelist_target))[-seq_len(nrow(edgelist_source))], ]

#Filter out clients active in the source dataset 
active_clients_source <-edgelist_source[edgelist_source$client %in% new_edges_ts$client,]

#Detect new edges between the prediction and target
new_edges_pt <-edgelist_pred[!duplicated(rbind(edgelist_target, edgelist_pred))[-seq_len(nrow(edgelist_target))], ]
#Filter out clients active in the source dataset
active_clients_target <-edgelist_target[edgelist_target$client %in% new_edges_pt$client,]


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


sample_a<-edgelist_source #after data frame
sample_b<-edgelist_target #before data frame
sample_c<-edgelist_pred #predict data frame

#sample_a<-sc_source #after data frame
#sample_b<-sc_target #before data frame
#sample_c<-sc_pred #predict data frame


g <- graph.data.frame(sample_a,directed = FALSE)
V(g)$type <- V(g)$name %in% sample_a[,1]
#is.bipartite(g) verify graph is bipartite
#get.edgelist(g, names=TRUE)# to verify this was done properly
adj_a <- get.adjacency(g, sparse = TRUE)


#isSymmetric(as.matrix(adj))
#calculate the eigenvalues and eigenvectors using ARPACK algorithm

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a %*% x) }
baev <- arpack(f2, sym=TRUE, options=list(n=vcount(g), nev=round(vcount(g)/40), ncv=round(vcount(g)/4+2),
                                          which="LM", maxiter=vcount(g)*12))
#plot(baev$values)
g_b <- graph.data.frame(sample_b,directed = FALSE)
V(g_b)$type <- V(g_b)$name %in% sample_b[,1]
#get.edgelist(g_b) #to verify this was done properly

adj_b <- get.adjacency(g_b, sparse = TRUE)
#isSymmetric(as.matrix(adj))
#calculate the eigenvalues and eigenvectors using ARPACK algorithm
f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b %*% x) }
baev_b <- arpack(f2, sym=TRUE, options=list(n=vcount(g_b), nev=round(vcount(g_b)/40), ncv=round(vcount(g_b)/4+2),
                                          which="LM", maxiter=vcount(g_b)*12))
u = Matrix(baev$vector,sparse = TRUE)
#ut = t(baev$vector) %*% adj_b %*% baev$vector
ut <- t(u) %*% adj_b 
ut <- ut %*% u
#image(ut)

f2 <- function(x, extra=NULL) { cat("."); as.vector(ut %*% x) }
baev_c <- arpack(f2, sym=FALSE, options=list(n=nrow(ut), nev=round(nrow(ut)/2), ncv=nrow(ut),
                                            which="LM", maxiter=nrow(ut)^2))

sort(baev_c$values,decreasing = FALSE)

g_c <- graph.data.frame(sample_c,directed = FALSE)
V(g_c)$type <- V(g_c)$name %in% sample_c[,1]
#get.edgelist(g_c) #to verify this was done properly
adj_c <- get.adjacency(g_c, sparse = TRUE)

df = as.data.frame(sort(as.numeric(baev$values[1:length(baev_c$values)]),decreasing = FALSE))
df$after = sort(as.numeric(baev_c$values[1:length(baev_c$values)]),decreasing = FALSE)
colnames(df) = c("source","target")
plot(df$source,df$target)


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

fit4 <- nls(formula = df$target ~ alpha*df$source, start=c(alpha=0.1),data = df)
summary(fit4)
sum(resid(fit4)^2)
confint(fit4)
new = data.frame(xdata = seq(min(df$source),max(df$source),len=length(df$source)))
lines(df$source,predict(fit3,newdata=new),col="brown")

#alpha <- 8.714e-11
#df$source <- alpha*df$source
#ggplot(data=df,aes(source,target))+geom_point()+
#  stat_smooth(method = "lm", formula = df$target ~ (I(df$source) + I(df$source^3)+I(df$source^5)+I(df$source^7)+I(df$source^9)),col = 'blue')+
#  stat_smooth(method = "lm", formula = df$target ~ df$source/(1-df$source^2),col = 'red')+
#  stat_smooth(method = "lm", formula = df$target ~ sinh(df$source),col='green')


#Now construct a prediction based on the best fit
flamda = matrix(NA,  ncol = length(baev_c$values),nrow = length(baev_c$values))
diag(flamda) = predict(fit1,baev_b$values)
flamda[is.na(flamda)] <- 0
flamda = Matrix(flamda,sparse = TRUE)
mp <- Matrix(baev_b$vectors[,1:length(baev_c$values)] %*% flamda %*% t(baev_b$vectors[,1:length(baev_c$values)]))
#set diagonal elements of prediction matrix to 0
dimnames(mp) <- dimnames(adj_c)
#image(mp)
#image(adj_c)

t1 <- as.vector(0)
for (i in 1: length(new_edges_pt$client))
{
  print(i)
  t1[i] <- mp[as.character(new_edges_pt$client[i]),as.character(new_edges_pt$merchant[i])]
}

mean(t1)
mean(t2)

hist(t1)

t2 <- as.vector(mp)
hist(t3)


#set all negative values to 0, all positive to 1
t<-(abs(mp)+mp)/(2*abs(mp))
diag(t) <- 0
image(t)


