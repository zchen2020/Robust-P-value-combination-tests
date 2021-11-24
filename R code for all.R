##################  Input ###############################
##########################################################
library("matrixStats")
#### matrix.type = "Expo" means AR(1) correlation structure and corresponds to Model 1 in Figure1.
#### matrix.type = "Poly" means polynomial decay and corresponds to Model 2 in Figure1.
#### matrix.type = "SiG"  means singular matrices and corresponds to Model 3 in Figure1.
matrix.type<-"SiG"

##########################################################
###################  Functions ###########################
##########################################################

test.all<-function(Pval){
  p.cau<-0.5-atan(mean(tan((0.5-Pval)*pi)))/pi # p-value from Cauchy
  p.min<-min(1,length(Pval)*min(Pval)) # p-value from MinP
  p.mcm<-min(1,2*min(p.cau,p.min)) # p-value from MCM 
  p.cmc<-0.5-atan(mean(c(tan((0.5-p.cau)*pi),tan((0.5-p.min)*pi))))/pi # p-value from CMC
  p.all<-c(p.cau, p.min, p.mcm,p.cmc)    
  return(p.all)
}



Sigma.chol.generate<-function(matrix.type,p,para){
  if (matrix.type=="Expo"){
    rho<-para
    Sigma<-matrix(0,nrow = p,ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        Sigma[i,j]<-rho^{abs(i-j)}
      }
    }
  }
  
  if (matrix.type=="Poly"){
    r<-para
    Sigma<-matrix(0,nrow = p,ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        if (i==j){
          Sigma[i,j]<-1
        }else{
          Sigma[i,j]<-1/(0.7+abs(i-j)^r)
        }
      }
    }
  }
  
  if (matrix.type=="SiG"){
    d<-para
    A<-matrix(NA,nrow=(p/5),ncol=p)
    for (i in 1:(p/5)){
      for (j in 1:p){
        A[i,j]<-1/(i-j+d)  
      }
    }
    A<-A%*%diag(1/sqrt(colSums(A^2)))
    Sigma.chol<-t(A) 
  }
  
  #### chol decomposition
  if (matrix.type!="SiG"){
    #suppressWarnings(Sigma.chol<-chol(Sigma,pivot = TRUE))
    Sigma.chol<-chol(Sigma,pivot = TRUE)
    Sigma.chol<-Sigma.chol[,order(attr(Sigma.chol, "pivot"))]
    Sigma.chol<-t(Sigma.chol)
  }
  
  return(Sigma.chol)
}

# simulate p-values and get the empirical power value:
pow.sim<-function(B=1e6,matrix.type="SiG",para=0.2,p=20,n.sig=4,n.conc=2,mu=0.6,alpha=0.01){
  Sigma.chol<-Sigma.chol.generate(matrix.type=matrix.type,p=p,para=para)
  pp<-dim(Sigma.chol)[2]
  p.tests.1l<-p.tests.1r<-p.tests.2<-matrix(NA,nrow=B,ncol=4)
  q<-rep(0,p)
  if (n.sig>0) {
    dir<-rep(1,n.sig) 
    sam.size<-rep(1,n.sig) 
    q<--mu*dir*sam.size 
  }
  
  for(i in 1:B) {
    X<-rnorm(pp)
    dim(X)<-c(pp,1)
    X<-Sigma.chol%*%X
    # sbuset<-sample(1:p, n.sig, replace=F)
    X[1:p]<-X[1:p]+q
    Pval.1l<-pnorm(X) # left-sided p-value
    Pval.1r<-1-Pval.1l # right-sided p-value
    Pval.2<-2*pnorm(-abs(X)) # 2-sided p-value
    p.tests.1l[i,]<-test.all(Pval.1l)
    p.tests.1r[i,]<-test.all(Pval.1r)
    p.tests.2[i,]<-test.all(Pval.2)
  } 
  pow.1l<-colMeans(p.tests.1l<alpha)
  pow.1r<-colMeans(p.tests.1r<alpha)
  pow.2<-colMeans(p.tests.2<alpha)
  pow<-rbind(pow.1l,pow.1r,pow.2)
  return(round(pow,2))
}


# example:

# Empirical Type I error rate:
pow.sim(B=1e6,matrix.type="Expo",para=0.5,p=10,n.sig=0,n.conc=0,mu=0,alpha=0.01)

# Empirical power:
pow.sim(B=1e4,matrix.type="Expo",para=0.5,p=10,n.sig=2,n.conc=1,mu=3,alpha=0.01)


############################################################################
# R code for real data application
############################################################################

# z-test for combining independent p-values:
z.test<-function(x){
  n<-length(x) 
  q<-sum(qnorm(x,lower.tail = FALSE))
  p<-pnorm(q/sqrt(n), lower.tail =  FALSE) 
  return(p)
}

# min-p test for combining independent p-values:
min.test<-function(x){
  n<-length(x) 
  p1=min(x)
  p<-1-(1-p1)^n
  return(p)
}

# Fisher test for combining independent p-values:
fisher.test<-function(x){
  n<-length(x) 
  q<-sum(qchisq(x,df=2,lower.tail = FALSE))
  p<-pchisq(q,df=2*n, lower.tail =  FALSE) 
  return(p)
}

# original data (OR and 95% CI)
a<-c(0.51,0.78,0.73,0.42,0.39,0.71,0.69,1.37,0.63,1.54,1.18,1.05)
b<-c(2.39,1.21,1.72,2.75,1.95,2.3,2.08,10.6,1.79,5.63,4.72,2.7)
yi<-log(c(1.11,0.97,1.13,1.08,0.88,1.28,1.19,3.82,1.06,2.95,2.36,1.68))
vi<-(log(b/a)/3.92)^2

## left 1-sided p-values:
p_1_l<-pnorm(yi/sqrt(vi))
# right 1-sided p-values:
p_1_r<-pnorm(-yi/sqrt(vi))

## 2-sided p-values:
p_2<-pchisq(yi^2/vi,df=1,lower.tail = FALSE)

## use z-test, minp-test and fisher-test to combine independent p-values from left, right one-sided tests and 2-sided test:
## combine independent p-values from left, right one-sided tests and 2-sided test:
p.1l.z<-z.test(p_1_l)
p.1r.z<-z.test(p_1_r)
p.2.z<-z.test(p_2)

p.1l.min<-min.test(p_1_l)
p.1r.min<-min.test(p_1_r)
p.2.min<-min.test(p_2)

p.1l.fisher<-fisher.test(p_1_l)
p.1r.fisher<-fisher.test(p_1_r)
p.2.fisher<-fisher.test(p_2)

## use CCT, MinP, MCM, and CMC to ombine the two combined p-values (using left and right one-sided p-values)
test.all(c(p.1l.z,p.1r.z))
test.all(c(p.1l.min,p.1r.min))
test.all(c(p.1l.fisher,p.1r.fisher))


