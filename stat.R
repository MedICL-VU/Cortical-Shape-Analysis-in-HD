### load dependencies
install.packages(c("Matrix"), repos="http://cran.us.r-project.org")
install.packages(c("multcomp"), repos="http://cran.us.r-project.org")
install.packages(c("R.matlab"), repos="http://cran.us.r-project.org")
install.packages(c("lme4"), repos="http://cran.us.r-project.org")
library(multcomp)
library(R.matlab)
library(lme4)

Sys.setenv('R_MAX_VSIZE'=32000000000)
Sys.getenv('R_MAX_VSIZE')

### get data
path <- "/Users/Zachary/Desktop/Vandy/Coursework/U_Third_Year/Summer_2020/CorticalShape_Zach"
#Y0 <- readMat(paste(path,"/data/y_sd.mat",sep=""))
#Y0 <- readMat(paste(path,"/data/y_ct.mat",sep=""))
Y0 <- readMat(paste(path,"/data/y_lgi.mat",sep=""))
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

### extract parameters for lmer
Duration <- demographics$age_at_scan
Sex <- demographics$gender
Subject <- demographics$external_id
Scanner <- demographics$scanner
for (sub in subject) {
  inds = subject == sub
  entry = min(Duration[inds])
  Duration[inds] <- Duration[inds] - entry
}
dx <- demographics$cap_e_group
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
models <- c()
for (i in 1:ncol(Y0)) {
  Y = Y0[,i]
  data = data.frame(Duration,Low,Med,High,Sex,Subject,Scanner,Y)
  M = lmer(Y ~ 1 + Duration + Low + Med + High + Sex + (1|Subject),data=data)
  models <- cbind(models,M)
  #write(summary(M),file="results_file",append=TRUE, sep="\t")
  if (i %% 10000 == 0) {
    sprintf("Iteration: %i", i)
  }
}

### conducting general linear hypothesis tests
htests <- c()
contr <- rbind("B0-B2" = c(1,0,-1,0,0,0,0), 
               "B0-B3" = c(1,0,0,-1,0,0,0),
               "B0-B4" = c(1,0,0,0,-1,0,0))
for (model in models) {
  test = glht(model,linfct=mcp(tension=contr))
  htests <- cbind(htests,test)
}



