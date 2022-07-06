#THE DATA CAN BE OBTAINED FROM:
# https://ag.purdue.edu/facai/data/
# https://www.gfbinitiative.org/data
#==========================================================================================#
#CODE FROM LIANG ET AL. 2016
#https://github.com/jingjingliang2018/GFB1
#==========================================================================================#
data <- subset(data, P>0)
data <- subset(data, S>0)

quantile(data$S,0.99996)
quantile(data$P,0.99996)

data1 <- subset(data,data$S<=270 & data$P<=533 & data$S >0 & data$P>0)

logP <- log(data1$P)
# jig coordinates to avoid duplicated values
Lon1 <- data1$Lon+ runif(length(data1$Lon),-0.0001,0.0001)
Lat1 <- data1$Lat+ runif(length(data1$Lat),-0.0001,0.0001)
data1 <- cbind.data.frame(data1, logP, Lat1, Lon1)
#==========================================================================================#
# MY CHANGES TO USE THE DATA FRAME
#==========================================================================================#
data1$logS<-log(data1$S)
df_RD<-list()
# The df_RD contains
# [[1]] <- GFBI data
# [[2]] <- parameteres names
# [[3]] <- formula
# ONLY TO COMAPRE WITH OLS
# [[4]] <- regression parameters obtained with the full data (OLS)
# [[5]] <- standard deviation of the regression parameters estimated with the full dataset(OLS)

df_RD[[1]]<-data1
df_RD[[2]]<-c("Intercept","log(S)","G","T3","C1","C3","PET","IAA","E")
df_RD[[3]]<-"logP~ logS + G + T3 + C1 + C3 + PET + IAA + E"

mod<-lm(as.formula(df_RD[[3]]), data = df_RD[[1]])

df_RD[[4]]<-coef(mod)
df_RD[[5]]<-sqrt(diag(vcov(mod)))

names(df_RD[[4]])<-df_RD[[2]]
names(df_RD[[5]])<-df_RD[[2]]

#save(df_RD , file = "empirical_data.RData")
