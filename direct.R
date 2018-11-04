## direct.R
## simulations with direct effect

rm(list=ls());
library(MASS);
source("fxns.R");

## **************************************************************
## simulation parameters
## **************************************************************
n <- 200;
pk <- 0.1; ## minor allele frequency
sig2 <- 4;
p <- 10; q <- 100; r <- 2; ## p=genes, q=snps, r = clinical covs

## transcript model: G~S+X relationship
set.seed(1);
Sigma <- cor(matrix(rnorm(p*10),nrow=10));
GamS <- matrix(sample(0:1,p*q,replace=TRUE,prob=c(0.5,0.5)),
               nrow=q,ncol=p); ## each S affects 3 genes on ave
GamS[GamS==1] <- sample(c(-1,1),sum(GamS==1),replace=TRUE)*
  runif(sum(GamS==1),0.05,1);
GamX <- matrix(sample(0:1,r*p,replace=TRUE,prob=c(0.5,0.5)),
                 nrow=r,ncol=p); ## each X affects 7 genes on ave
GamX[GamX==1] <- sample(c(-1,1),sum(GamX==1),replace=TRUE)*
  runif(sum(GamX==1),0.05,1);

## outcome model with direct effect
## 1. linear in genes and clinical covs
## 2. correct gene set
## 3. direct effect
## 4. no measurement error
aG <- sample(c(-1,1),p,replace=TRUE)*runif(p,0.05,0.7);
aX <- rnorm(r);
## to be least favorable to integrative method, make direct
## effect such that standard gwas effect sizes all get increased
aS <- sign(drop(GamS%*%aG))*0.75; aS <- c(aS,0);

## add row to GamS, corresponding to unimportant SNP, by making
## the row orthogonal to aG
gamS <- rnorm(p);
gamS <- gamS-sum(gamS*aG)/sum(aG^2)*aG;
GamS <- rbind(GamS,gamS);

## **************************************************************
## simulate linear model and case control
## **************************************************************
## store results, .l = linear, .c = case control
## dimensions: results (beta, variance, z, p), snp, integrative
## or standard or overlap
res.l <- res.c <- array(NA,dim=c(4,q+1,3));

## run simulation for a given random seed
env <- 1;
set.seed(env);

## ==============================================================
## linear model
## ==============================================================
S <- matrix(sample(0:2,n*(q+1),replace=TRUE,
                   prob=c((1-pk)^2,2*pk*(1-pk),pk^2)),
            nrow=n);
X <- matrix(rnorm(n*r),ncol=r);
G <- S%*%GamS+X%*%GamX+mvrnorm(n,rep(0,p),sig2*Sigma);
Y <- G%*%aG+S%*%aS+X%*%aX+rnorm(n,0,2);

## integrative
res.l[,,1] <- apply(S,2,function(s)
                  {
                      fit.i <- igwas(Y,G,s,X,family="gaussian",
                                     direct = TRUE);
                      return(unlist(fit.i));
                  });

## standard gwas (gives t instead of z values);
res.l[,,2] <- apply(S,2,function(s)
                    {
                        fit.s <- summary(lm(Y~s+X))$coefficients;
                        return(c(fit.s[2,1],fit.s[2,2]^2,fit.s[2,3],fit.s[2,4]));
                    });

## overlap
genes.Y <- which(p.adjust(apply(G,2,function(g)
                                { summary(lm(Y~g+X))$coefficients[2,4]; }),
                          method="fdr")<=0.05);
## only return p-values
if(length(genes.Y)==0){ res.l[4,,3] <- rep(1,q+1); } else
{
    res.l[4,,3] <- apply(S,2,function(s)
        {
            ps.s <- apply(G,2,function(g){ summary(lm(g~s+X))$coefficients[2,4]; })
            genes.s <- which(p.adjust(ps.s,method="fdr")<=0.05);
            if(length(genes.s)==0){ return(1); } else
            {
                o <- length(intersect(genes.s,genes.Y));
                return(1-phyper(o-1,length(genes.Y),p-length(genes.Y),length(genes.s)));
            }
        });
}

## ==============================================================
## case control
## ==============================================================
## generate lots of data then sample cases and controls
nn <- 5000;
S <- matrix(sample(0:2,nn*(q+1),replace=TRUE,
                   prob=c((1-pk)^2,2*pk*(1-pk),pk^2)),
            nrow=nn);
X <- matrix(rnorm(nn*r),ncol=r);
G <- S%*%GamS+X%*%GamX+mvrnorm(nn,rep(0,p),sig2*Sigma);
Y <- sapply(expit(-5.8+G%*%aG+S%*%aS+X%*%aX),function(x)
            { sample(0:1,1,prob=c(1-x,x)); });
mean(Y);
## marginal prevalence of ~31%
P <- 0.31;

## sample cases and controls
keep <- c(which(Y==1)[1:(n/2)],which(Y==0)[1:(n/2)]);
S <- S[keep,];
X <- X[keep,];
G <- G[keep,];
Y <- Y[keep];

## integrative
res.c[,,1] <- apply(S,2,function(s)
                  {
                      fit.i <- igwas.cc(Y,P,G,s,X,family="binomial",
                                        direct = TRUE);
                      return(unlist(fit.i));
                  });

## standard gwas
res.c[,,2] <- apply(S,2,function(s)
                    {
                        fit.s <- summary(glm(Y~s+X,family="binomial"))$coefficients;
                        return(c(fit.s[2,1],fit.s[2,2]^2,fit.s[2,3],fit.s[2,4]));
                    });

## overlap
genes.Y <- which(p.adjust(apply(G,2,function(g)
                                { summary(glm(Y~g+X,family="binomial"))$coefficients[2,4]; }),
                          method="fdr")<=0.05);
## only return p-values
if(length(genes.Y)==0){ res.c[4,,3] <- rep(1,q+1); } else
{
    ## need secondary phenotype methods to calculate p-values
    w1 <- P*n/sum(Y==1); w0 <- (1-P)*n/sum(Y==0);
    ws <- Y; ws[Y==1] <- w1; ws[Y==0] <- w0;
    res.c[4,,3] <- apply(S,2,function(s)
        {
            ps.s <- apply(G,2,function(g)
                { summary(lm(g~s+X,weights=ws))$coefficients[2,4]; })
            genes.s <- which(p.adjust(ps.s,method="fdr")<=0.05);
            if(length(genes.s)==0){ return(1); } else
                {
                    o <- length(intersect(genes.s,genes.Y));
                    return(1-phyper(o-1,length(genes.Y),p-length(genes.Y),length(genes.s)));
                }
        });
}

## **************************************************************
## save results
## **************************************************************
rm(list=setdiff(ls(),c("res.l","res.c","env")));
save.image(paste("results/direct",env,".RData",sep=""));
