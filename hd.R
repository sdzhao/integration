## hd.R
## simulation for easy case: correctly specified, low dim

rm(list=ls());
library(MASS);
source("fxns.R");

## **************************************************************
## simulation parameters
## **************************************************************
n <- 400;
pk <- 0.1; ## minor allele frequency
sig2 <- 16;
q <- 10000; p <- 2*q; r <- 2; ## p=genes, q=snps, r = clinical covs

## transcript model: G~S+X relationship
set.seed(1);
## each snp has two cis genes, X's also play a role
gamS <- matrix(rnorm(2*q,0,5),ncol=2);
GamX <- matrix(sample(0:1,r*p,replace=TRUE,prob=c(0.5,0.5)),
                 nrow=r,ncol=p); ## each X affects 7 genes on ave
GamX[GamX==1] <- sample(c(-1,1),sum(GamX==1),replace=TRUE)*
  runif(sum(GamX==1),0.05,1);
## match snps (1st col) with genes (2nd col)
cispairs <- cbind(kronecker(1:q,rep(1,2)),1:p);

## high-dimensional outcome model
## 1. linear in genes and clinical covs
## 2. high-dimensional gene set
## 3. no direct effect
## 4. no measurement error
aG <- sample(0:1,p,replace=TRUE,prob=c(0.999,0.001));
aG[aG==1] <- sample(c(-1,1),sum(aG),replace=TRUE)*
   runif(sum(aG),1,5);
aX <- rnorm(r);

## **************************************************************
## simulate linear model and case control
## **************************************************************
## store p-values, .l = linear, .c = case control
res.l.i <- res.c.i <- rep(NA,nrow(cispairs)); ## integrative
res.l.g <- res.c.g <- rep(NA,q); ## standard gwas
##res.l.o <- res.c.o <- rep(NA,q); ## overlap

## run simulation for a given random seed
env <- 1;
set.seed(env);

## ==============================================================
## linear model
## ==============================================================
S <- matrix(sample(0:2,n*q,replace=TRUE,
                   prob=c((1-pk)^2,2*pk*(1-pk),pk^2)),
            nrow=n);
X <- matrix(rnorm(n*r),ncol=r);
G <- matrix(NA,nrow=n,ncol=p);
for(i in 1:q)
{ G[,(2*(i-1)+1):(2*i)] <- matrix(S[,i],ncol=1)%*%gamS[1,]; }
G <- G+X%*%GamX+
    t(replicate(n,arima.sim(list(order=c(1,0,0),ar=0.5),p)))*sqrt(sig2)*sqrt(1-0.5^2);
Y <- G%*%aG+X%*%aX+rnorm(n,0,2);

## integrative
res.l.i <- apply(cispairs,1,function(x)
                 { igwas(Y,G[,x[2]],S[,x[1]],X,family="gaussian")$p; });

## standard gwas
res.l.g <- apply(S,2,function(s){ summary(lm(Y~s+X))$coefficients[2,4]; });

## these take too long
## ## overlap
## genes.Y <- which(p.adjust(apply(G,2,function(g)
##                                 { summary(lm(Y~g+X))$coefficients[2,4]; }),
##                           method="fdr")<=0.05);
## ## only return p-values
## if(length(genes.Y)==0){ res.l.o <- rep(1,q); } else
## {
##     res.l.o <- apply(S,2,function(s)
##         {
##             ps.s <- apply(G,2,function(g){ summary(lm(g~s+X))$coefficients[2,4]; })
##             genes.s <- which(p.adjust(ps.s,method="fdr")<=0.05);
##             if(length(genes.s)==0){ return(1); } else
##             {
##                 o <- length(intersect(genes.s,genes.Y));
##                 return(1-phyper(o-1,length(genes.Y),p-length(genes.Y),length(genes.s)));
##             }
##         });
## }

## ==============================================================
## case control
## ==============================================================
## generate lots of data then sample cases and controls
nn <- 2000;
S <- matrix(sample(0:2,nn*q,replace=TRUE,
                   prob=c((1-pk)^2,2*pk*(1-pk),pk^2)),
            nrow=nn);
X <- matrix(rnorm(nn*r),ncol=r);
G <- matrix(NA,nrow=nn,ncol=p);
for(i in 1:q)
{ G[,(2*(i-1)+1):(2*i)] <- matrix(S[,i],ncol=1)%*%gamS[1,]; }
G <- G+X%*%GamX+
    t(replicate(nn,arima.sim(list(order=c(1,0,0),ar=0.5),p)))*sqrt(sig2)*sqrt(1-0.5^2);
Y <- sapply(expit(-16+G%*%aG+X%*%aX),function(x)
            { sample(0:1,1,prob=c(1-x,x)); });
## marginal prevalence of ~31%
P <- 0.31;

## sample cases and controls
keep <- c(which(Y==1)[1:(n/2)],which(Y==0)[1:(n/2)]);
S <- S[keep,];
X <- X[keep,];
G <- G[keep,];
Y <- Y[keep];

## integrative
res.c.i <- apply(cispairs,1,function(x)
                 { igwas.cc(Y,P,G[,x[2]],S[,x[1]],X,family="binomial")$p; });

## standard gwas
res.c.g <- apply(S,2,function(s)
                 { summary(glm(Y~s+X,family="binomial"))$coefficients[2,4]; });

## these take too long
## ## overlap
## genes.Y <- which(p.adjust(apply(G,2,function(g)
##                                 { summary(glm(Y~g+X,family="binomial"))$coefficients[2,4]; }),
##                           method="fdr")<=0.05);
## ## only return p-values
## if(length(genes.Y)==0){ res.c.o <- rep(1,q); } else
## {
##     res.c.o <- apply(S,2,function(s)
##         {
##             ps.s <- apply(G,2,function(g)
##                           { summary(glm(g~s+X,family="binomial"))$coefficients[2,4]; })
##             genes.s <- which(p.adjust(ps.s,method="fdr")<=0.05);
##             if(length(genes.s)==0){ return(1); } else
##             {
##                 o <- length(intersect(genes.s,genes.Y));
##                 return(1-phyper(o-1,length(genes.Y),p-length(genes.Y),length(genes.s)));
##             }
##         });
## }

## **************************************************************
## save results
## **************************************************************
rm(list=setdiff(ls(),c("res.l.i","res.l.g","res.c.i","res.c.g","env")));
save.image(paste("results/hd",env,".RData",sep=""));
