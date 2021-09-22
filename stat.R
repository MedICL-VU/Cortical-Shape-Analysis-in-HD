### load dependencies
install.packages(c("Matrix"), repos="http://cran.us.r-project.org")
install.packages(c("multcomp"), repos="http://cran.us.r-project.org")
install.packages(c("R.matlab"), repos="http://cran.us.r-project.org")
install.packages(c("lme4"), repos="http://cran.us.r-project.org")
install.packages(c("foreach"), repos="http://cran.us.r-project.org")
install.packages(c("doParallel"), repos="http://cran.us.r-project.org")
library(Matrix)
library(multcomp)
library(R.matlab)
library(lme4)
library(foreach)
library(doParallel)

Sys.setenv('R_MAX_VSIZE'=32000000000)
Sys.getenv('R_MAX_VSIZE')

### get data
cmeasure = "lgi" # {ct_smooth=6,sd,lgi}
path <- getwd()
Y0 <- readMat(sprintf("%s/data/y_%s.mat",path,cmeasure))
Y0 <- Y0[1]$Y0
demographics <- read.csv(paste(path,"/data/demographics.csv",sep=""))
env <- readMat(paste(path,"/env/environment.mat",sep=""))

### remove unknown cap_e_group values
valid <- demographics$cap_e_group==0 | demographics$cap_e_group==1 | demographics$cap_e_group==2 | demographics$cap_e_group==3
demographics <- demographics[valid,]
Y0 <- Y0[valid,]

### get the subjects that we want
subj <- (demographics$cag > 0 & demographics$cag < 37) | demographics$cag > 37
Y0 <- Y0[subj,]
demographics <- demographics[subj,]

# correct SD correlation with ticv
if (cmeasure == "sd") {
  Y0 <- diag(1/demographics$ticv) %*% Y0
}

### extract parameters for lmer
Duration <- demographics$age_at_scan
Age <- demographics$age_at_scan
Sex <- demographics$gender
Subject <- demographics$external_id
Scanner <- demographics$scanner
for (sub in Subject) {
  inds = Subject == sub
  entry = min(Duration[inds])
  Duration[inds] <- Duration[inds] - entry
}
dx <- demographics$cap_e_group
Cntrl <- dx
Cntrl[dx == 0] = 1
Cntrl[dx != 0] = 0
Low <- dx
Low[dx == 1] = 1
Low[dx != 1] = 0
Med <- dx
Med[dx == 2] = 1
Med[dx != 2] = 0
High <- dx
High[dx == 3] = 1
High[dx != 3] = 0

### fitting the linear mixed model
Y0[,sum(abs(Y0)) == 0] <- runif(length(Y0[,sum(abs(Y0))==0]))*.Machine$double.eps # prevent numerical instability

####################################################################################################
### DEBUG
####################################################################################################
Y = Y0[,1]
data = data.frame(Duration,Age,Sex,Low,Med,High,Subject,Y)
M = lmer(Y ~ 0 + Duration + Age + Sex + Cntrl + Low + Med + High + (1 | Subject),data=data)

#write(summary(M),file="results_file",append=TRUE, sep="\t")

### conducting general linear hypothesis tests
contr <- rbind("B3-B4" = c(0,0,0,1,-1,0,0), 
               "B3-B5" = c(0,0,0,1,0,-1,0),
               "B3-B6" = c(0,0,0,1,0,0,-1))

#test = summary(glht(M,linfct=mcp(tension=contr)))
test = summary(glht(M,linfct=contr,df=3), test=Chisqtest())
chi = test[["test"]][["SSH"]]
coefs <- test[["model"]]@beta
  
####################################################################################################

### MAIN LOOP
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

chis_betas <- foreach(i=1:ncol(Y0), .combine=rbind, .packages=c('lme4','multcomp')) %dopar% {
  Y = Y0[,i]
  data = data.frame(Duration,Age,Sex,Low,Med,High,Subject,Y)
  M = lmer(Y ~ 0 + Duration + Age + Sex + Cntrl + Low + Med + High + (1 | Subject),data=data)
  
  #write(summary(M),file="results_file",append=TRUE, sep="\t")
  
  ### conducting general linear hypothesis tests
  contr <- rbind("B3-B4" = c(0,0,0,1,-1,0,0), 
                 "B3-B5" = c(0,0,0,1,0,-1,0),
                 "B3-B6" = c(0,0,0,1,0,0,-1))
  
  #test = summary(glht(M,linfct=mcp(tension=contr)))
  test = summary(glht(M,linfct=contr,df=3), test=Chisqtest())
  chi = test[["test"]][["SSH"]]
  coefs <- test[["model"]]@beta
  c(chi,coefs)
}
#stop cluster
stopCluster(cl)

chis <- chis_betas[,1]
betas <- chis_betas[,2:ncol(chis_betas)]
writeMat(sprintf("%s/data/chis_betas_%s.mat",path,cmeasure),chis=chis,betas=betas)


