require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)
require(SnowballC)
require(lsa)
require(scales)
require(reshape2)
require(corpcor)
options(java.parameters = "-Xmx4g")

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

queryC <- paste("select * from pp_oap_sing_es_t.transactions_clients")
client_list <- dbGetQuery(conn, queryC)

queryM <- paste("select * from pp_oap_sing_es_t.transactions_merchants")
merchant_list <- dbGetQuery(conn, queryM)

print("closing db connection")
dbDisconnect(conn)


t1<-edgelist_source[,2:3]
t2<-edgelist_target[,2:3]
t3<-edgelist_pred[,2:3]


#Define number of eigenvalues to be used with similarity metrics
q <- 20

find_eigenvalues <- function(df)
{
  g <- graph.data.frame(df,directed = FALSE)
  adj <- get.adjacency(g, sparse = TRUE,type=c("both"))
  f2 <- function(x, extra=NULL) { cat("."); as.vector(adj %*% x) }
  baev <- arpack(f2, sym=TRUE, options=list(n=vcount(g), nev=q, ncv=q+3,
                                            which="LM", maxiter=vcount(g)*12))
  
  return (list(baev,ecount(g)))
}

tlist <- list(t1,t2,t3)

#Find edge count and list of eigenvalues 
#se denotes spectral evolution
#se is an object that contains eigenvalue decomposition for 
#t1-x timepoints together with the respective edge count
se <- lapply(tlist,find_eigenvalues)

##Spectral evolution of the network
eigen_values <- Matrix(nrow=length(se), ncol = q)
edge_count <- list()

for (i in 1 : length(se))
{  
  for (j in 1: q)
  {
    eigen_values[i,j] <- se[[i]][[1]]$values[j]
  }
  edge_count[i] <- se[[i]][[2]]
}

plot(edge_count,eigen_values[,1],xlab = 'Edge count', ylab = 'Eigenvalues',type = 'b',col = hcl(h = 15, c = 100, l = 85),ylim = c(min(eigen_values),max(eigen_values))
     ,main = 'Spectral evolution. Dominant values closer to red',xlim=c(as.numeric(edge_count[1]),as.numeric(edge_count[length(edge_count)])),pch=20,cex = 1)

for (i in 2:q)
{
  lines(edge_count,eigen_values[,i],xlab = 'Edge count', ylab = 'Eigenvalues',type = 'b',col = hcl(h = 0+15*i, c = 100, l = 85),pch=20,cex = 1)
}

#Eigenvector evolution
#Computes chart of eigenvector evolution
sim_k_l <- Matrix(nrow=length(se), ncol = q)
edge_count <- list()

for (i in 1 : length(se))
{  
  for (j in 1: q)
  {
    sim_k_l[i,j] <- abs(cosine(as.vector(se[[i]][[1]]$vectors[,j]),as.vector(se[[1]][[1]]$vectors[,j])))
  }
  edge_count[i] <- se[[i]][[2]]
}

plot(edge_count,sim_k_l[,1],xlab = 'Edge count', ylab = 'Similarity (sim(k,k))',type = 'b',col = hcl(h = 15, c = 100, l = 85),ylim=c(0,1)
     ,main = 'Eigenvector evolution Dominant vectors closer to red',xlim=c(as.numeric(edge_count[1]),as.numeric(edge_count[length(edge_count)])),pch=20,cex = 1)

for (i in 2:q)
{
  lines(edge_count,sim_k_l[,i],xlab = 'Edge count', ylab = 'Similarity (sim(k,k))',type = 'b',col = hcl(h = 0+15*i, c = 100, l = 85),pch=20,cex = 1)
}

#Eigenvector stability
#Computes stability of all the eigenvectors from time to time T
#substitute with the required X for t and Y for T in 
#se[[X]][[1]] and se[[Y]][[1]]
sim_k_k <- Matrix(nrow=q, ncol = q)

for (i in 1:q)
{
  for (j in 1:q)
  {
    sim_k_k[i,j] <- abs(cosine(as.vector(se[[1]][[1]]$vectors[,i]),as.vector(se[[2]][[1]]$vectors[,j])))
  }
}

image(sim_k_k)

#Spectral diagonality

sd <- function(df1,df2)
{
  #compute matrix of new edges - i.e. B
  g_1 <- graph.data.frame(df1,directed = FALSE)
  adj_1 <- get.adjacency(g_1, sparse = TRUE,type=c("both"))
  row.order <- rownames(adj_1)
  col.order <- colnames(adj_1)
  g_2 <- graph.data.frame(df2,directed = FALSE)
  adj_2 <- get.adjacency(g_2, sparse = TRUE,type=c("both"))
  adj_2 <- adj_2[row.order,col.order]
  ne <- adj_2 - adj_1
  ne <- ne[row.order,col.order]
  f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_1 %*% x) }
  baev_1 <- arpack(f2, sym=TRUE, options=list(n=vcount(g_1), nev=q, ncv=q+3,
                                              which="LM", maxiter=vcount(g_1)*12))
  #compute multiplication of eigenvalues of A and matrix B
  b_eigen_1 <- t(baev_1$vectors) %*% ne %*% baev_1$vectors
  #hist(as.numeric(b_eigen))
  image(b_eigen_1)
}

sd(t1,t2)


sample_a <- t1
sample_b <- t2
sample_c <- t3

#set the number of eigenvalues to be used in the prediction
r=20

g_a <- graph.data.frame(sample_a,directed = FALSE)
adj_a <- get.adjacency(g_a, sparse = TRUE,type=c("both"))

row.order <- rownames(adj_a)
col.order <- colnames(adj_a)

g_b <- graph.data.frame(sample_b,directed = FALSE)
adj_b <- get.adjacency(g_b, sparse = TRUE,type=c("both"))
adj_b <- adj_b[row.order,col.order]


#Check the development of the prediction eigenvalues in comparison to target
g_c <- graph.data.frame(sample_c,directed = FALSE)
adj_c <- get.adjacency(g_c, sparse = TRUE,type=c("both"))
adj_c <- adj_c[row.order,col.order]


f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a %*% x) }
baev_a <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b %*% x) }
baev_b <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                            which="LM", maxiter=vcount(g_a)*12))


#compute matrix of new edges - i.e. B
ne <- adj_b - adj_a
hist(as.numeric(ne))
ne <- ne[row.order,col.order]
#ne[ne >= 1] <- 1
#ne[ne < 0] <- 0

#normalise the spectra
#nlamda <- baev_a$values * abs(baev_b$values[1])/abs(baev_a$values[1])

#without normalisation
#nlamda <- baev_a$values


#compute multiplication of eigenvalues of A and matrix B
b_eigen <- t(baev_a$vectors) %*% ne %*% baev_a$vectors
#hist(as.numeric(b_eigen))
image(b_eigen)

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
flamda_a_polyn = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_polyn) = predict(a_polyn,new.dfb)
mp_a_polyn <- baev_b$vectors %*% flamda_a_polyn %*% t(baev_b$vectors)
dimnames(mp_a_polyn) <- dimnames(adj_a)
mp_a_polyn <- Matrix(mp_a_polyn,sparse=TRUE)

flamda_a_exp = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_exp) = predict(a_exp,new.dfb)
mp_a_exp <- baev_b$vectors %*% flamda_a_exp %*% t(baev_b$vectors)
dimnames(mp_a_exp) <- dimnames(adj_a)
mp_a_exp <- Matrix(mp_a_exp,sparse=TRUE)

flamda_a_neu = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b$values))
diag(flamda_a_neu) = predict(a_neu,new.dfb)
mp_a_neu <- baev_b$vectors %*% flamda_a_neu %*% t(baev_b$vectors)
dimnames(mp_a_neu) <- dimnames(adj_a)
mp_a_neu <- Matrix(mp_a_neu,sparse=TRUE)

#Compute the test matrix
netest <- adj_c - adj_b
hist(as.numeric(netest))
#image(netest) #new edges
#substitute any multiple edges with just 1 edge
netest[netest >= 1] <- 1
netest[netest < 0] <- 0


#Rank reduction
adj_b_rr <- Matrix(nrow=nrow(adj_b),ncol=ncol(adj_b),sparse=TRUE)
lambda = Matrix(0,  ncol = r,nrow = r,sparse=TRUE)
diag(lambda) = baev_b$values
adj_b_rr <- baev_b$vectors %*% lambda %*% t(baev_b$vectors)

#calculate jaccard similarity measure
adj_b_jacc <- Matrix(similarity.jaccard(g_b),sparse=TRUE)

#calculate adamic adar similarity measure
adj_b_adar <- Matrix(similarity.invlogweighted(g_b),sparse=TRUE)

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


#######################################Normalised adjacency matrix

#normalise matrix adj_a
adj_a_d <- degree(adj_a)

#take squre root of the diagonal of the degree matrix
adj_a_ds <- Matrix(data = 0,nrow=nrow(adj_a),ncol=nrow(adj_a))

diag(adj_a_ds)<-1/sqrt(diag(adj_a_d))

#compute normalised adj_a matrix
adj_a_n <- adj_a_ds%*%adj_a%*%adj_a_ds

#normalise matrix adj_b
adj_b_d <- degree(adj_b)

adj_b_ds <- Matrix(data = 0,nrow=nrow(adj_a),ncol=nrow(adj_a))

diag(adj_b_ds)<-1/sqrt(diag(adj_b_d))

#compute normalised adj_b matrix
adj_b_n <- adj_b_ds%*%adj_b%*%adj_b_ds

#t <- (adj_a_n + t(adj_a_n))

#ne_n <- adj_b + t(adj_b) - (adj_a + t(adj_a))

ne_n <- adj_b_n - adj_a_n


f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a_n %*% x) }
baev_a_n <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+10,
                                              which="LM", maxiter=vcount(g_a)*20))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b_n%*% x) }
baev_b_n <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+10,
                                              which="LM", maxiter=vcount(g_a)*20))


#Spectral evolution does not apply to normalised matrix N
#compute multiplication of eigenvalues of A and matrix B
b_eigen_n <- t(baev_a_n$vectors) %*% ne_n %*% baev_a_n$vectors
#hist(as.numeric(b_eigen))
image(b_eigen_n)

df = as.data.frame(baev_a_n$values) #the independent variable
df$flamda = as.numeric(diag(b_eigen_n)) #the dependent function
colnames(df) = c("source","target")
plot(x=df$source,y=df$target)


#Normalised path counting via polynomial
n_polyn <-nls(formula = target ~ I(alpha*source) + I(beta*source^2) + I(gamma*source^3) + I(delta*source^4) + I(delta*source^5),start = c(alpha=0.01,beta=0.01,gamma=0.01,delta=0.01), data = df)
summary(n_poly)
sum(resid(n_poly)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_poly,new)
points(x,y,col = "green")

#Normalised exponential kernel 
n_exp <- nls(formula = target ~ I(alpha*source)+I(((alpha*source)^2)/2)+I(((alpha*source)^3)/6+I(((alpha*source)^4)/24)),start=c(alpha=0.01),data=df)
summary(n_exp)
sum(resid(n_exp)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_exp,new)
points(x,y,col = "red")

#Normalised Neumann kernel
n_neu <- nls(formula = target ~ I(alpha*source)+I((alpha*source)^2)+I((alpha*source)^3)+I((alpha*source)^4)+I((alpha*source)^5),start=c(alpha=0.01),data = df)
summary(n_neu)
sum(resid(n_neu)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(n_neu,new)
points(x,y,col = "blue")

#Now construct a prediction based on the best fit
flamda_n_polyn = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_n$values))
diag(flamda_n_polyn) = predict(n_polyn,new.dfb)
mp_n_polyn <- baev_b_n$vectors %*% flamda_n_polyn %*% t(baev_b_n$vectors)
dimnames(mp_n_polyn) <- dimnames(adj_a)
mp_n_polyn <- Matrix(mp_n_polyn,sparse = TRUE)

flamda_n_exp = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_n$values))
diag(flamda_n_exp) = predict(n_exp,new.dfb)
mp_n_exp <- baev_b_n$vectors %*% flamda_n_exp %*% t(baev_b_n$vectors)
dimnames(mp_n_exp) <- dimnames(adj_a)
mp_n_exp <- Matrix(mp_n_exp,sparse = TRUE)


flamda_n_neu = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_n$values))
diag(flamda_n_neu) = predict(n_neu,new.dfb)
mp_n_neu <- baev_b_n$vectors %*% flamda_n_neu %*% t(baev_b_n$vectors)
dimnames(mp_n_neu) <- dimnames(adj_a)
mp_n_neu <- Matrix(mp_n_neu,sparse = TRUE)


###########################################Laplacian matrix
#compute Laplacian
adj_a_l <- adj_a_d - adj_a
adj_b_l <- adj_b_d - adj_b

#normalise the laplacian
adj_a_l_n <- adj_a_ds%*%adj_a_l%*%adj_a_ds
adj_b_l_n <- adj_b_ds%*%adj_b_l%*%adj_b_ds

#non-normalised Laplacian kernel
adj_b_l_plus <- Matrix(pseudoinverse(adj_b_l) ,sparse = TRUE)


#normalised commute-time kernel (normalised Laplacian kernel)
adj_b_l_n_plus <- Matrix(pseudoinverse(adj_b_l_n),sparse = TRUE)

#new edges Laplacian matrix
ne_l <- adj_b_l - adj_a_l

#compute smallest eigenvalues for Laplacian here
f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_a_l %*% x) }
baev_a_l <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                              which="LM", maxiter=vcount(g_a)*40))

f2 <- function(x, extra=NULL) { cat("."); as.vector(adj_b_l %*% x) }
baev_b_l <- arpack(f2, sym=TRUE, options=list(n=vcount(g_a), nev=r, ncv=r+3,
                                              which="LM", maxiter=vcount(g_a)*40))


b_eigen_l <- t(baev_a_l$vectors) %*% ne_l %*% baev_a_l$vectors
#hist(as.numeric(b_eigen))
image(b_eigen_l)

df = as.data.frame(baev_a_l$values) #the independent variable
df$flamda = as.numeric(diag(b_eigen_l)) #the dependent function
colnames(df) = c("source","target")
plot(x=df$source,y=df$target)

#Regularised Laplacian kernel
l_reg <-nls(formula = target ~ 1/(1+alpha*I(source)),start = c(alpha=0.00001), data = df)
summary(l_reg)
sum(resid(l_reg)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(l_reg,new)
points(x,y,col = "green")

#Heat diffusion kernel
l_hd <- nls(formula = target ~ I(-alpha*source)+I(((-alpha*source)^2)/2)+I(((-alpha*source)^3)/6+I(((-alpha*source)^4)/24)),start=c(alpha=0.01),data=df)
summary(l_hd)
sum(resid(l_hd)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(l_hd,new)
points(x,y,col = "red")

#Generalised Laplacian kernel
l_gen <-nls(formula = target ~ 1+1/(alpha*(1-I(source)))+1/(alpha*(1-I(source))^2),start = c(alpha=0.1,beta=0.1), data = df)
summary(l_gen)
sum(resid(l_gen)^2)
new = data.frame(source = df$source)
x=df$source
y=predict(l_gen,new)
points(x,y,col = "blue")

#Now construct a prediction based on the best fit
flamda_l_reg = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_l$values))
diag(flamda_l_reg) = predict(l_reg,new.dfb)
mp_l_reg <- baev_b_l$vectors %*% flamda_l_reg %*% t(baev_b_l$vectors)
dimnames(mp_l_reg) <- dimnames(adj_a)
mp_l_reg <- Matrix(mp_l_reg,sparse = TRUE)

flamda_l_hd = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_l$values))
diag(flamda_l_hd) = predict(l_hd,new.dfb)
mp_l_hd <- baev_b_l$vectors %*% flamda_l_hd %*% t(baev_b_l$vectors)
dimnames(mp_l_hd) <- dimnames(adj_a)
mp_l_hd <- Matrix(mp_l_hd,sparse = TRUE)

flamda_l_gen = matrix(0,  ncol = r,nrow = r)
new.dfb <- data.frame(source = as.numeric(baev_b_l$values))
diag(flamda_l_gen) = predict(l_gen,new.dfb)
mp_l_gen <- baev_b_l$vectors %*% flamda_l_gen %*% t(baev_b_l$vectors)
dimnames(mp_l_gen) <- dimnames(adj_a)
mp_l_gen <- Matrix(mp_l_gen,sparse = TRUE)

######################################### EVALUATION

####Calculate prediction quality
pred_a_polyn <- prediction(as.vector(mp_a_polyn),as.vector(netest))
pred_a_exp <- prediction(as.vector(mp_a_exp),as.vector(netest))
pred_a_neu <- prediction(as.vector(mp_a_neu),as.vector(netest))
pred_rr <- prediction(as.vector(adj_b_rr),as.vector(netest))

pred_n_exp <- prediction(as.vector(mp_n_exp),as.vector(netest))
pred_n_neu <- prediction(as.vector(mp_n_neu),as.vector(netest))

pred_jacc <- prediction(as.vector(adj_b_jacc),as.vector(netest))
pred_adar <- prediction(as.vector(adj_b_adar),as.vector(netest))
pred_ppa <- prediction(as.vector(adj_b_ppa),as.vector(netest))

pred_b_plus <- prediction(as.vector(adj_b_l_plus),as.vector(netest))
pred_b_l_n_plus <- prediction(as.vector(adj_b_l_n_plus),as.vector(netest))

pred_l_reg <- prediction(as.vector(mp_l_reg),as.vector(netest))
pred_l_hd <- prediction(as.vector(mp_l_hd),as.vector(netest))

######################################ROC curves

###ADJACENCY MATRIX ROC CURVES
perf_a_polyn <- performance(pred_a_polyn, "tpr", "fpr")
par(lty=1)
plot(perf_a_polyn,col = "black")

perf_a_exp <- performance(pred_a_exp, "tpr", "fpr")
par(new = TRUE)
plot(perf_a_exp,col = "light grey")

perf_a_neu <- performance(pred_a_neu, "tpr", "fpr")
par(new=TRUE)
plot(perf_a_neu,col = "grey")

perf_rr <- performance(pred_rr,"tpr","fpr")
par(new=TRUE)
plot(perf_rr,col = "dark grey")

###NORMALISED ADJACENCY MATRIX ROC CURVES
perf_n_exp <- performance(pred_n_exp, "tpr", "fpr")
par(new = TRUE,lty = 2)
plot(perf_n_exp,col = "orange")

perf_n_neu <- performance(pred_n_neu, "tpr", "fpr")
par(new = TRUE)
plot(perf_n_exp,col = "yellow")

###LOCAL METRICS ROC CURVES
perf_jacc <- performance(pred_jacc,"tpr","fpr")
par(new=TRUE,lty = 3)
plot(perf_jacc, col = "light blue")

perf_adar <- performance(pred_adar,"tpr","fpr")
par(new=TRUE)
plot(perf_adar,col = "blue")

perf_ppa <- performance(pred_ppa,"tpr","fpr")
par(new=TRUE)
plot(perf_ppa, col = "dark blue")


###LAPLACIAN METRICS ROC CURVES
perf_b_plus <- performance(pred_b_plus,"tpr","fpr")
par(new=TRUE,lty = 4)
plot(perf_b_plus,col = "red")

perf_b_l_n_plus <- performance(pred_b_l_n_plus,"tpr","fpr")
par(new=TRUE)
plot(perf_b_l_n_plus,col = "darkred")

perf_l_reg <- performance(pred_l_reg,"tpr","fpr")
par(new=TRUE)
plot(perf_l_reg,col = "red2")

perf_l_hd <- performance(pred_l_hd,"tpr","fpr")
par(new=TRUE)
plot(perf_l_hd,col = "orangered")



############################Precision recall curves

precision_recall_a_polyn <- performance(pred_a_polyn, "prec", "rec")
plot(precision_recall_a_polyn,ylim = c(0,1),col = "black")

precision_recall_a_exp <- performance(pred_a_exp, "prec", "rec")
par(new=TRUE)
plot(precision_recall_a_exp,ylim = c(0,1),col = "yellow")

precision_recall_a_neu <- performance(pred_a_neu, "prec", "rec")
par(new=TRUE)
plot(precision_recall_a_neu,ylim=c(0,1),col="pink")

precision_recall_rr <- performance(pred_rr, "prec", "rec")
par(new=TRUE)
plot(precision_recall_rr,ylim=c(0,1),col="brown")

precision_recall_exp <- performance(pred_n_exp, "prec", "rec")
par(new = TRUE)
plot(precision_recall_exp,ylim=c(0,1),col = "orange")

precision_recall_n_neu <- performance(pred_n_neu, "prec", "rec")
par(new = TRUE)
plot(precision_recall_n_neu,ylim=c(0,1),col = "light green")

precision_recall_jacc <- performance(pred_jacc, "prec", "rec")
par(new=TRUE)
plot(precision_recall_jacc,ylim=c(0,1),col="blue")

precision_recall_adar <- performance(pred_adar, "prec", "rec")
par(new=TRUE)
plot(precision_recall_adar,ylim=c(0,1),col="green")

precision_recall_ppa <- performance(pred_ppa, "prec", "rec")
par(new=TRUE)
plot(precision_recall_ppa,ylim = c(0,1),col = "red")

precision_recall_perf_b_plus <- performance(pred_b_plus, "prec", "rec")
par(new=TRUE)
plot(precision_recall_perf_b_plus,ylim = c(0,1),col = "cyan")

precision_recall_b_l_n_plus <- performance(pred_b_l_n_plus,"prec","rec")
par(new=TRUE)
plot(precision_recall_b_l_n_plus,ylim = c(0,1),col = "purple")

precision_recall_l_reg <- performance(pred_l_reg,"prec","rec")
par(new=TRUE)
plot(precision_recall_l_reg,ylim = c(0,1),col = "dark red")

precision_recall_l_hd <- performance(pred_l_hd,"prec","rec")
par(new=TRUE)
plot(precision_recall_l_hd,ylim = c(0,1),col = "dark green")



############################Phi values

precision_f_b_l_n_plus <- performance(pred_b_l_n_plus,"f")

precision_f_a_polyn <- performance(pred_a_polyn,"f")

precision_prbe_b_l_n_plus <- performance(pred_b_l_n_plus,"prbe")


##########################################Extract predictions

cutoffs <- data.frame(cut=perf_rr@alpha.values[[1]], fpr=perf_rr@x.values[[1]], 
                      tpr=perf_rr@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
head(cutoffs)
head(subset(cutoffs, fpr < 0.10))

cutoffs <- data.frame(cut=precision_recall_rr@alpha.values[[1]], recall=precision_recall_rr@x.values[[1]], 
                      precision=precision_recall_rr@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$precision, decreasing=TRUE),]
head(cutoffs)
head(subset(cutoffs, recall > 0.5))

###Spectral extrapolation

pred <- adj_b_rr
dimnames(pred) <- dimnames(adj_a)

pred[pred>=1] <- 1
pred[pred<1] <- 0

ind <- which(pred == 1)

for (i in 1:length(ind))
{
  k <- arrayInd(ind[i],dim(pred))
  rownames(pred)[k[,1]]
  colnames(pred)[k[,2]]
  print(paste('Predicted negative score',pred[ind[i]],' for the the following link'))
  print(mapply(`[[`, dimnames(pred), k))
  print(paste('In reality the following value was observed in the adjacency matrix'))
  print(adj_c[rownames(pred)[k[,1]],colnames(pred)[k[,2]]])
  print(paste('New edges matrix had the following value'))
  print(netest[rownames(pred)[k[,1]],colnames(pred)[k[,2]]])
  print("")
}



#Find the largest connected component of a graph
g_a <- graph.data.frame(edgelist_source[,2:3],directed = FALSE)
gclust = clusters(g_a)
x <- which.max(sizes(gclust))
g_a1 <- induced.subgraph(g_a, which(membership(gclust) == x))

g_b <- graph.data.frame(edgelist_target[,2:3],directed = FALSE)
g_b1 <- induced.subgraph(g_b, which(membership(gclust) == x))

g_c <- graph.data.frame(edgelist_pred[,2:3],directed = FALSE)
g_c1 <- induced.subgraph(g_c, which(membership(gclust) == x))

vcount(g_a1)
vcount(g_b1)
vcount(g_c1)

ecount(g_a1)
ecount(g_b1)
ecount(g_c1)

edgelist_source <- get.data.frame(g_a)
edgelist_target <- get.data.frame(g_b)
edgelist_pred <- get.data.frame(g_c)

t1 <- edgelist_source
t2 <- edgelist_target
t3 <- edgelist_pred


