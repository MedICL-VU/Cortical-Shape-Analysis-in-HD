library(multcomp)
hsb2 <- read.csv("https://stats.idre.ucla.edu/stat/data/hsb2.csv")

m1 <- lm(read ~ socst + factor(ses) * factor(female), data = hsb2)
summary(m1)