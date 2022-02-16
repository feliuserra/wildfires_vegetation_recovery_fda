#####################
# Load libraries
rm(list=ls())

library(refund)
library(fda.usc)
# library(dplyr)

library(mgcv)

# Set working directory and load data
#####################

setwd('~/BSC/luminosity/fda_wildfires/')
#setwd('../')
#load('.RData')
df_trends <- read.csv('data/processed/FDA_analysis/wildfires_recoveries_trend_with_5pre_7post.csv')
#df_trends_ife <- read.csv('data/processed/wildfires_recoveries_trend_with_5pre_7post_with_IFE.csv')
df_trends$X <- NULL
#df_trends_ife$X <- NULL
#df_trends_wout_extrapolate <- read.csv('data/processed/wildfires_recoveries_trend_with_5pre_7post_removing_extrapolation.csv')
dict_treatments <- jsonlite::fromJSON('data/processed/FDA_analysis/dict_wildfires_times.json')
df_covariates <- read.csv('data/processed/gdf_covariates_final.csv')
df_covariates$X <- NULL
df_covariates$avg_altitude <- df_covariates$avg_slope
df_covariates$avg_slope <- NULL
df_covariates$logAcres <- log(df_covariates$Acres)

dim(df_trends)
#dim(df_trends_ife)
#dim(df_trends_wout_extrapolate)
length(dict_treatments)
dim(df_covariates)


####### load(file='SOFR_session_13042021.RData')

####################

# Preprocess data
#####################

df_wildfires <- df_trends
df_wildfires_covariates <- df_covariates[df_covariates$ID_new %in% paste(names(df_trends),'_A', sep=''),]

dim(df_wildfires)
dim(df_wildfires_covariates)

sum(paste(names(df_wildfires), '_A', sep='') == as.character(df_wildfires_covariates$ID_new))

list_landcovers <- list()
count <- 1
for (x in df_wildfires_covariates$landcover){
  if (df_wildfires_covariates$landcover[count] == 'shrub_scrub'){
    list_landcovers[count] <- 'shrub_scrub'
  }
  else if (df_wildfires_covariates$landcover[count] == 'evergreen_forest'){
    list_landcovers[count] <- 'evergreen_forest'
  }
  else if (df_wildfires_covariates$landcover[count] == 'grassland_herbaceous'){
    list_landcovers[count] <- 'grassland_herbaceous'
  }
  else {
    list_landcovers[count] <- 'other_landcover'
  }
  count <- count + 1
}

df_wildfires_covariates$landcover_grouped <- as.factor(unlist(list_landcovers))
  
#####################

# Load data as fdata and make basic exploratory data analysis
#####################
# fda.usc library
years_argvals <- seq(from=0, to=7, length.out=182)
my.recovery <- fdata(mdata=t(df_wildfires), 
                     argvals=years_argvals,
                     names=list(main='LSR curves',
                                xlab='Years since wildfire',
                                ylab=expression(paste(delta, ' NDVI'))))


plot(my.recovery, col='grey', xlim=c(0,7))
a1<-func.mean(my.recovery)
lines(a1, lty=1, lwd=2, col=2)
abline(h=0, lty=2)

###########
df_lm <- df_wildfires_covariates
df_lm$ID <- NULL
df_lm$Fire_ID <- NULL

df_lm_orig <- df_wildfires_covariates
df_lm_orig$ID <- NULL
df_lm_orig$Fire_ID <- NULL

ind <- sapply(df_lm, is.numeric)
df_lm[,ind] <- apply(df_lm[,ind], 2, scale)
summary(df_lm)

df_lm_orig$Fire_Type <- relevel(as.factor(df_lm_orig$Fire_Type), ref = "WF")

df_lm$Fire_Type <- relevel(as.factor(df_lm$Fire_Type), ref = "WF")


# Make FOSR regression
####################
df_lm_orig$response <- my.recovery$data
df_lm$response <- my.recovery$data


################
# Fitting all the univariate regression in loop
#
list_regres <- list( c("lat","lon"),  "avg_altitude", 
                     "Year", "StartMonth",
  "logAcres", c("landcover_grouped","entropy_la"),
  c("avg_ndvi_5_year_before", "std_ndvi_5_year_before"), 
  "bi", "tmmx", "pr", "srad")
n.regs <- length(list_regres)

Pctg.Dev.Expl <- matrix(nrow = n.regs, ncol=5)
models.aic <- matrix(nrow = n.regs, ncol=5)

list_models <- vector(n.regs,mode="list")
for (i in 1:n.regs) list_models[[i]] <- vector(5, mode="list")

sm.par <- rep("",n.regs)
sm.par[4]<-", bs='cc', k=7"

# c()
for (i in 1:n.regs){
  if (i==6){# c("landcover_grouped","entropy_la")
    regr <- list_regres[i]
    formu.char <- paste0("response ~ c(",regr[[1]][1],") + c(",regr[[1]][2],")")
    formu <- as.formula(formu.char)
  }else{
    if (i %in% c(1,7)){
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(",regr[[1]][1],") + c(",regr[[1]][2],")")
      formu <- as.formula(formu.char)
    }else{
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(",regr,")")
      formu <- as.formula(formu.char)
    }
  }
  print(paste(i, " ", formu.char))
  sofr_i <- pffr(formu, data=df_lm, yind=years_argvals)
  Pctg.Dev.Expl[i,1] <- summary(sofr_i)$dev.expl*100
  models.aic[i,1] <- sofr_i$aic
  list_models[[i]][[1]] <- sofr_i
}

# c(s())
for (i in 1:n.regs){
  if (i==6){# c("landcover_grouped","entropy_la")
    regr <- list_regres[i]
    formu.char <- paste0("response ~ c(",regr[[1]][1],") + c(s(",regr[[1]][2],"))")
    formu <- as.formula(formu.char)
  }else{
    if (i %in% c(1,7)){
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(s(",regr[[1]][1],")) + c(s(",regr[[1]][2],"))")
      formu <- as.formula(formu.char)
    }else{
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(s(",regr,sm.par[i],"))")
      formu <- as.formula(formu.char)
    }
  }
  print(paste(i, " ", formu.char))
  sofr_i <- pffr(formu, data=df_lm, yind=years_argvals)
  Pctg.Dev.Expl[i,2] <- summary(sofr_i)$dev.expl*100
  models.aic[i,2] <- sofr_i$aic
  list_models[[i]][[2]] <- sofr_i
}

# no c(), no s()
for (i in 1:n.regs){
  if (i==6){# c("landcover_grouped","entropy_la")
    regr <- list_regres[i]
    formu.char <- paste0("response ~ ",regr[[1]][1]," + ",regr[[1]][2])
    formu <- as.formula(formu.char)
  }else{
    if (i %in% c(1,7)){
      regr <- list_regres[i]
      formu.char <- paste0("response ~ ",regr[[1]][1]," + ",regr[[1]][2])
      formu <- as.formula(formu.char)
    }else{
      regr <- list_regres[i]
      formu.char <- paste0("response ~ ",regr)
      formu <- as.formula(formu.char)
    }
  }
  print(paste(i, " ", formu.char))
  sofr_i <- pffr(formu, data=df_lm, yind=years_argvals)
  Pctg.Dev.Expl[i,3] <- summary(sofr_i)$dev.expl*100
  models.aic[i,3] <- sofr_i$aic
  list_models[[i]][[3]] <- sofr_i
}

# c(s(regressor)) + regressor 
for (i in 1:n.regs){
  if (i==6){# c("landcover_grouped","entropy_la")
    regr <- list_regres[i]
    formu.char <- paste0("response ~ c(",regr[[1]][1],") + c(s(",regr[[1]][2],")) + ",regr[[1]][2])
    formu <- as.formula(formu.char)
  }else{
    if (i %in% c(1,7)){
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(s(",regr[[1]][1],")) + ",
                           regr[[1]][1], " + c(s(",regr[[1]][2],")) + ",
                           regr[[1]][2])
      formu <- as.formula(formu.char)
    }else{
      regr <- list_regres[i]
      formu.char <- paste0("response ~ c(s(",regr,sm.par[i],")) + ",regr)
      formu <- as.formula(formu.char)
    }
  }
  print(paste(i, " ", formu.char))
  sofr_i <- pffr(formu, data=df_lm, yind=years_argvals)
  Pctg.Dev.Expl[i,4] <- summary(sofr_i)$dev.expl*100
  models.aic[i,4] <- sofr_i$aic
  list_models[[i]][[4]] <- sofr_i
}

# s()
for (i in 1:n.regs){
  if (i==6){# c("landcover_grouped","entropy_la")
    regr <- list_regres[i]
    formu.char <- paste0("response ~ ",regr[[1]][1]," + s(",regr[[1]][2],")")
    formu <- as.formula(formu.char)
  }else{
    if (i %in% c(1,7)){
      regr <- list_regres[i]
      formu.char <- paste0("response ~ s(",regr[[1]][1],") + s(",regr[[1]][2],")")
      formu <- as.formula(formu.char)
    }else{
      regr <- list_regres[i]
      formu.char <- paste0("response ~ s(",regr,sm.par[i],")")
      formu <- as.formula(formu.char)
    }
  }
  print(paste(i, " ", formu.char))
  sofr_i <- pffr(formu, data=df_lm, yind=years_argvals)
  Pctg.Dev.Expl[i,5] <- summary(sofr_i)$dev.expl*100
  models.aic[i,5] <- sofr_i$aic
  list_models[[i]][[5]] <- sofr_i
}

################

apply(Pctg.Dev.Expl,1,which.max)
#  [1] 4 5 4 4 4 4 5 4 5 4 5

apply(models.aic,1,which.min)
#  [1] 4 5 4 4 4 4 5 4 5 4 5

# Conclusion: 
# The selected univariate models are the same using "Pctg.Dev.Expl" or "models.aic"

bmatrix = function(x, digits=NULL, ...) {
  library(xtable)
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE)
  passed_args = list(...)
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  cat("\\begin{bmatrix}\n",
      do.call(print.xtable, calling_args),
      "\\end{bmatrix}\n")
}


bmatrix(Pctg.Dev.Expl, digits=2)

################

# Prueba para ver si funciona anova()
i<- 9 
regr <- list_regres[i]

formu.char <- paste0("response ~ c(",regr,")")
formu <- as.formula(formu.char)
print(paste(i, " ", formu.char))
model.9.1 <- pffr(formu, data=df_lm, yind=years_argvals)

formu.char <- paste0("response ~ ",regr)
formu <- as.formula(formu.char)
print(paste(i, " ", formu.char))
model.9.3 <- pffr(formu, data=df_lm, yind=years_argvals)

formu.char <- paste0("response ~ c(s(",regr,")) + ",regr)
formu <- as.formula(formu.char)
print(paste(i, " ", formu.char))
model.9.4 <- pffr(formu, data=df_lm, yind=years_argvals)

anova(model.9.1,model.9.3)
# anova(model.9.1,model.9.3,test="F") # da error

# Conclusión: anova.gam() no hace F-tests con las salidas de pffr()

#############################
# 25 de abril (antes de usar el AIC)
# Elección de los términos individuales que entrarán en el modelo global
#
# Regla: Para cada regresor, elegimos la forma en que entra (c(), c(s()), etc.)
# que mayor "Pctg. Dev. Expl." ten ga, 
# salvo que la inmediatamente anterior (y más simple) esté 
# a menos de un 1% de "Pctg. Dev. Expl."
# Ese término más simple entra en el modelo a menos que 
# el inmediatamente anterior esté a menos de un 1% de "Pctg. Dev. Expl.",
# y así sucesivamente.
#
# Con este criterio, el modelo elegido lleva:
# - Todas las variables de la forma c(s(variable))
# - salvo las variables 7 y 10, que entran como c(s(variable)) + variable

# A full model
sofr_full<- pffr(response ~ 
                    c(s(lat)) + c(s(lon)) +
                    c(s(avg_altitude)) + 
                    c(s(Year)) + 
                    c(s(StartMonth, bs='cc', k=7)) +
                    c(s(logAcres)) + 
                    c(landcover_grouped) + c(s(entropy_la)) +
                    avg_ndvi_5_year_before + std_ndvi_5_year_before + 
                    c(s(avg_ndvi_5_year_before)) + c(s(std_ndvi_5_year_before)) + 
                    c(s(bi)) + 
                    c(s(tmmx)) +
                    c(s(pr)) + pr + 
                    c(s(srad)),
                  data=df_lm,
                  yind = years_argvals)
summary(sofr_full) # Deviance explained =   72.9%
plot(sofr_full, pages=1, scale=0)

plot(sofr_full, pages=4, scale=-1)


coefs <- coef(sofr_full)

dataframe_constant_values <- as.data.frame(coefs$pterms)

dataframe_values <- data.frame(1:100)
dataframe_values_s <- data.frame(1:8000)
dataframe_values_c_s <- data.frame(1:100)
for (i in names(coefs$smterms)){
  if (startsWith(i, 'c(s(')){
    eval(parse(text=paste("dataframe_values_c_s$`", i, "`<- coefs$smterms$`", i, "`$value", sep='')))
    eval(parse(text=paste("dataframe_values$`", i, "_se`<- coefs$smterms$`", i, "`$se", sep='')))
    eval(parse(text=paste("dataframe_values$`index", i, "` <- coefs$smterms$`", i, "`$x", sep='')))
  } else {
    eval(parse(text=paste("dataframe_values$`", i, "`<- coefs$smterms$`", i, "`$value", sep='')))
    eval(parse(text=paste("dataframe_values$`", i, "_se`<- coefs$smterms$`", i, "`$se", sep='')))
    eval(parse(text=paste("dataframe_values$`index", i, "` <- coefs$smterms$`", i, "`$x", sep='')))
  }
}

list_dataframes <- c(dataframe_constant_values, dataframe_values, dataframe_values_c_s)
write.csv(list_dataframes, 'results/covariates_effects_processed_for_plot.csv')

save.image('Final_Object_pffrs_27042021.RData')

#############################