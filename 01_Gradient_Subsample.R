# LIST OF FUNCTIONS
# GS_SD_01, GS_GD_01 and GS_RD_01 are the Gradient Subsampling functions 
# that generate the pseudo-parameters.
# 
# GS_02 calculates the sample mean and sample standard deviation
# GS_03 estimate the parameters and the standard deviation for the full size

#==========================================================================================#
# DATA

#==========================================================================================#
#NECESSARY LIBRARIES
library(tidyverse)
library(lme4)
library(nlme)
#==========================================================================================#
#STEP 1: PSEUDO-PARAMETERS
#==========================================================================================#
# INPUT
# full_df: df_LinR, df_LogR of df_NLinR
# b_subset: numerical vector with the subsample sizes
# repetitions: number of pseudo-replicates 
# parameters: vector with the names of the regression parameters
# formula: formula to use in the regression analysis (see formulas for each regression model)
# regression.model: type of regression model, the options are 
# "LinR"(lm) -> "y~x1+x2"
# "LogR" (glm, family = binomial(link = "logit")) -> "y~x1+x2"
# "NLinR" (nls, start = list(a = 1, b = 1)) -> "y~(a*x)/(b+x)"
# OUTPUT: 
# OUT_1[[1]]<- pseudo-parameters 
# OUT_1[[2]]<- self reported parameter SE
#==========================================================================================#
GS_SD_01<-function(full_df, repetitions, b_subset, parameters, formula, regression.model){
  
  f.01<-as.formula(formula)
  
  sub_samples<- list(sub_param = list(), 
                     sub_SE = list())
  
  sub_samples_est<- list()
  sub_samples_SE<-list()
  
  N<-nrow(full_df)
  
  for (j in 1:length(b_subset)){
    
    b_subsample<- b_subset[j]
    sub_samples_df<- data.frame()
    
    
    sub_samples_est[[j]]<- data.frame(matrix(ncol = length(parameters)))
    colnames(sub_samples_est[[j]])<- parameters
    
    sub_samples_SE[[j]]<- data.frame(matrix(ncol = length(parameters)))
    colnames(sub_samples_SE[[j]])<- parameters
    
    
    for(i in 1:repetitions){
      #DRAW THE SAMPLES AND MAKE A DF 
      sub_samples_df<- full_df[sample(nrow(full_df), size = b_subsample, replace = F), ]
      
      #DO A REGRESSION 
      
      if(regression.model == "LinR"){
        sub_samples_mod<-lm(f.01, data = sub_samples_df) 
      }
      
      if(regression.model == "LogR"){
        sub_samples_mod<-lm(f.01, data = sub_samples_df, family = binomial(link = "logit")) 
      }
      
      if(regression.model == "NLinR"){
        sub_samples_mod<-nls(f.01, start = list(a = 1, b = 1),data = sub_samples_df) 
      }
      
      
      #SAVE THE COEFCIENTS
      sub_samples_est[[j]][i,]<- coef(sub_samples_mod)
      
      #SAVE THE SE OF EACH ESTIMATED PARAMETER
      sub_samples_SE[[j]][i,]<- sqrt(diag(vcov(sub_samples_mod)))
      
    }
    
    names(sub_samples_est)[j]<- b_subsample
    names(sub_samples_SE)[j]<- b_subsample
    
    
  }
  
  #Regression parameters
  sub_samples[[1]] <- sub_samples_est
  
  #SE of the Regression parameter
  sub_samples[[2]] <- sub_samples_SE
  
  return(sub_samples)
}

#==========================================================================================#
# INPUT
# full_df: LMM1, LMM2 or LMM3
# b_subset: numerical vector with the subsample sizes
# repetitions: number of pseudo-replicates 
# parameters: vector with the names of the regression parameters
# group.level: group factor in the data
# formula: formula to use in the regression analysis (specify the random effects according to the package lme4)
# use the following formula according to the LMM data
# LMM1 -> "y~x1+x2+(1|grp)"
# LMM2 -> "y~x1+x2+(1|grp)"
# LMM3 -> "y~x1+x2+(0+x2|grp)"
# OUTPUT: 
# OUT_1[[1]]<- pseudo-parameters 
# OUT_1[[2]]<- random-effects pseudo parameters
# OUT_1[[3]]<- self reported parameter SE
#==========================================================================================#


GS_GD_01<-function(full_df, repetitions, b_subset, parameters, group.levels ,formula){
  
  
  f.01<- as.formula(formula)
  
  sub_samples<- list(sub_param = list(),
                     sub_ranef = list(),
                     sub_SE = list())
  
  sub_samples_est<- list()
  sub_samples_SE<-list()
  sub_samples_ranef<-list()
  
  
  
  for (j in 1:length(b_subset)){
    
    b_subsample<- b_subset[j]
    sub_samples_df<- data.frame()
    
    # #MODEL
    
    #1 ESTIMATED PARAMETERS
    sub_samples_est[[j]]<- data.frame(matrix(ncol = length(parameters)))
    colnames(sub_samples_est[[j]])<- parameters
    
    #2 ESTIMATED PARAMETERES SE
    sub_samples_SE[[j]]<- data.frame(matrix(ncol = length(parameters)))
    colnames(sub_samples_SE[[j]])<- parameters
    
    #3 RANDOM EFFECTS
    sub_samples_ranef[[j]]<-data.frame(matrix(ncol = length(group.levels)))
    colnames(sub_samples_ranef[[j]])<-group.levels
    
    
    for(i in 1:repetitions){
      
      
      #DRAW THE SAMPLES AND MAKE A DF 
      sub_samples_df<- full_df[sample(nrow(full_df), size = b_subsample, replace = F), ]
      
      #LMM WITH THE MADE SUBSET
      sub_samples_mod <-lmer(f.01 , data = sub_samples_df , REML = T)
      
      
      #1 SAVE ESTIMATED PARAMETER
      sub_samples_est[[j]][i,]<-fixef(sub_samples_mod)
      
      #2 SAVE SE OF THE ESTIMATED PARAMETER
      sub_samples_SE[[j]][i,]<-sqrt(diag(vcov(sub_samples_mod)))
      
      #3 SAVE THE RANDOM EFFECTS
      #IF NOT ALL EFFECTS ARE PRESENT, I HAVE TO KEEP THE POSITIONS OF THE ONES THAT ARE PRESENT
      re.new<-rep(0,length(group.levels))
      
      grp.lev<-as.data.frame(ranef(sub_samples_mod))[ ,"grp"]
      re<-as.data.frame(ranef(sub_samples_mod))[ ,"condval"]
      
      missing<-which(group.levels%in% grp.lev == F) 
      present<-which(group.levels%in% grp.lev == T)
      
      re.new[present]<-re
      re.new[missing]<-NA
      
      sub_samples_ranef[[j]][i,]<-re.new
      
      
      
    }
  names(sub_samples_est)[j]<-b_subset[j]  
  names(sub_samples_SE)[j]<-b_subset[j]
  names(sub_samples_ranef)[j]<-b_subset[j] 
    
  }
  
  
  sub_samples[[1]]<-sub_samples_est
  sub_samples[[2]]<-sub_samples_ranef
  sub_samples[[3]]<-sub_samples_SE
  
  return(sub_samples)
}


#==========================================================================================#
# INPUT
# full_df: df_RD[[1]]
# repetitions: number of pseudo-replicates 
# b_subset: numerical vector with the subsample sizes
# parameters: vector with the names of the regression parameters
# formula: formula to use in the regression analysis
# regression.model: type of regression model, the options are 
# "OLS"(lm)
# "GLS"(gls)
# OUTPUT: 
# OUT_1[[1]]<- pseudo-parameters 
#==========================================================================================#

GS_RD_01<-function(full_df, repetitions, b_subset, parameters, formula, regression.model){
  
  f.01<-as.formula(formula)
  
  sub_samples<- list(sub_param = list())
  
  sub_samples_est<- list()
  
  N<-nrow(full_df)
  
  for (j in 1:length(b_subset)){
    
    b_subsample<- b_subset[j]
    sub_samples_df<- data.frame()
    
    
    sub_samples_est[[j]]<- data.frame(matrix(ncol = length(parameters)))
    colnames(sub_samples_est[[j]])<- parameters
    
    
    
    for(i in 1:repetitions){
      #DRAW THE SAMPLES AND MAKE A DF 
      
      sub_samples_df<- full_df[sample(nrow(full_df), size = b_subsample, replace = F), ]
      
      #DO A REGRESSION 
      
      if(regression.model == "OLS"){
        sub_samples_mod<-lm(f.01, data = sub_samples_df) 
      }
      
      if(regression.model == "GLS"){
        tryCatch({
        sub_samples_mod<-gls(f.01 , data = sub_samples_df , method="ML", 
                             corr = corSpher(form = ~ Lon + Lat, nugget = TRUE), 
                             control = glsControl(singular.ok=TRUE)) 
      },error=function(e){})
      }
      
      
      #SAVE THE COEFCIENTS
      sub_samples_est[[j]][i,]<- coef(sub_samples_mod)
      

      
    }
    
    names(sub_samples_est)[j]<- b_subsample
    
    
  }
  
  #Regression parameters
  sub_samples[[1]] <- sub_samples_est
  
  
  return(sub_samples)
}

#==========================================================================================#
# EXAMPLES: THE DATAFRAMES SHOLD BE LOADED FIRST

#==========================================================================================#
# EXAMPLE SIMPLE DATA
#==========================================================================================#
R<-100
b<-seq(10,100, by = 10)
parameters.names<-letters[1:3]
df<-df_LinR[[1]]

OUT_1<-GS_SD_01(full_df = df, repetitions = R, b_subset = b,
                parameters = parameters.names, formula = "y~x1+x2", regression.model = "LinR")


#==========================================================================================#
# EXAMPLE GROUPED DATA
#==========================================================================================#
R<-100
b<-seq(10,100, by = 10)
parameters.names<-letters[1:3]
df<-df_LMM1[[1]]
grp<-names(table(df$grp))

OUT_1<-GS_GD_01(full_df = df, repetitions = R, b_subset = b,
                parameters = parameters.names, group.levels = grp ,formula = "y~x1+x2+(1|grp)")

#==========================================================================================#
# EXAMPLE REAL DATA
#==========================================================================================#
R<-100
b<-seq(10,100, by = 10)
parameters.names<-df_RD[[2]]
df<-df_RD[[1]]

grp<-names(table(df$grp))

OUT_1<-GS_RD_01(full_df = df, repetitions = R, b_subset = b,
                parameters = parameters.names, regression.model = "OLS",
                formula = df_RD[[3]])

#==========================================================================================#
#STEP 2: MEAN AND STANDARD DEVIATION

# INPUT
# OUT_1[[1]]: pseudo-parameters
# parameters: vector with the names of the regression parameters
# OUTPUT: 
# OUT_2 <- data frame with the statistic mean and standard deviation
#==========================================================================================#

GS_02<-function(out.1, parameters){
  
  pseudo_parameters<-do.call(rbind, out.1)
  pseudo_parameters$size<-trunc(as.numeric(rownames(pseudo_parameters)))
  row.names(pseudo_parameters)<-NULL
  
  df_1<-pseudo_parameters%>%
    gather(as.factor(parameters), key = "parameter", value = "value")%>%
    group_by(size, parameter)%>%
    summarise(stat.mean = mean(value), stat.sd = sd(value))
  
  
  return(as.data.frame(df_1))
}

OUT_2<-GS_02(out.1 = OUT_1[[1]], parameters = parameters.names)


#==========================================================================================#
#STEP 3: GS estimates 

# INPUT
# OUT_2:data frame with the statistic mean and standard deviation
# parameters: vector with the names of the regression parameters
# N: length of the full data set
# OUTPUT: 
# OUT_3 <- data frame with the statistic estimates 
#==========================================================================================#


GS_03<-function(out.2, parameters, N){
  
  C.b<-rep(0, length(parameters))
  
  for (i in 1:length(parameters)){
    mod<-lm(log(stat.sd)~ offset(-0.5*log(size)), data = out.2[out.2$parameter == parameters[i],])
    C.b[i]<-exp(coef(mod))
  }
  
  stat.sd.b<-C.b/sqrt(N)
  
  df_1<-out.2%>%
    group_by(parameter)%>%
    summarise(stat.value.b = weighted.mean(stat.mean, size))
  
  df_2<-cbind(df_1, stat.sd.b, C.b)
  
  return(df_2)
}

OUT_3<-GS_03(out.2 = OUT_2, parameters = parameters.names , N = nrow(df))

df_RD[[5]]
#==========================================================================================#
#STANDARD DEVIATION RANDOM EFFECT

# INPUT (ONLY FROM GROUPED DATA)
# OUT_1[[2]]: pseudo-pseudo estimates random-effects
# group.level: vector with the ID of the grouping factor
# OUTPUT: 
# OUT_re <- data frame with the random-effect sample standard deviation calculated for each subsample size
#==========================================================================================#

GS_RE<-function(out.1, group.levels){
  
  
  
  
  pseudo_parameters<-do.call(rbind, out.1)
  pseudo_parameters$size<-trunc(as.numeric(rownames(pseudo_parameters)))
  row.names(pseudo_parameters)<-NULL
  
  df_1<-pseudo_parameters%>%
    gather(as.factor(group.levels), key = "group", value = "value")%>%
    group_by(size, group)%>%
    summarise(stat.mean = mean(value, na.rm = T))
  
  re.sd<-df_1%>%
    group_by(size)%>%
    summarise(re.sd = sqrt((1/(length(group.levels)-1))*sum((stat.mean)^2)) ) 
  
  
  return(as.data.frame(re.sd))
}

OUT_re<-GS_RE(out.1 = OUT_1[[2]], group.levels = grp)


