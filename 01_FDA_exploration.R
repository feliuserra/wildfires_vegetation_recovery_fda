
library(refund)
library(fda)
library(fda.usc)

setwd('~/BSC/luminosity/fda_wildfires/')
#df <- read.csv('data/processed/FDA_analysis/wildfires_effects_with_5pre_10post.csv')
df <- read.csv('../google_earth_engine_lab/data/processed/FDA_analysis/wildfires_effects_trend_with_5pre_10post.csv')
df_trends <- read.csv('data/processed/FDA_analysis/wildfires_effects_trend_with_7pre_7post.csv')
df_trends <- read.csv('../google_earth_engine_lab/data/processed/FDA_analysis/wildfires_effects_trend_with_5pre_10post.csv')
dict_treatments <- jsonlite::fromJSON('data/processed/FDA_analysis/dict_wildfires_times.json')
df_covariates <- read.csv('../../luminosity/google_earth_engine_lab/data/interim/gdf_now_landcover.csv')

dim(df)
dim(df_trends)

list_w <- list()
list_w_covariates <- list()

list_years <- list()
list_month <- list()
list_day <- list()
list_fire_type <- list()
list_acres <- list()
list_state <- list()
list_landcover <- list()
list_landcover_entropy <- list()

for (wildfire in names(df_trends)[2:length(df_trends)]){
  list_w[[wildfire]] <- df_trends[dict_treatments[[wildfire]]:
              (dict_treatments[[wildfire]]+26*10),
            wildfire]
  list_w_covariates[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, c('Year', 'StartMonth', 
                                                                                      'StartDay', 'Fire_Type', 'Acres', 
                                                                                      'State', 'landcover', 'entropy_la')]
  list_years[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'Year']
  list_month[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'StartMonth']
  list_day[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'StartDay']
  list_fire_type[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'Fire_Type']
  list_acres[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'Acres']
  list_state[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'State']
  list_landcover[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'landcover']
  list_landcover_entropy[[wildfire]] <- df_covariates[df_covariates$Fire_ID == wildfire, 'entropy_la']
}

df_wildfires <- data.frame(list_w)
order(names(df_wildfires))
df_wildfires <- df_wildfires[, order(names(df_wildfires))]
df_wildfires_covariates <- data.frame(list_w_covariates)


years_argvals <- seq(from=0, to=261 / (365 / 14), length.out=261)
my.recovery <- fdata(mdata=t(df_wildfires), 
                     argvals=years_argvals,
                     names=list(main='LSR curves',
                                xlab='Years since wildfire',
                                ylab='NDVI'))

plot.fdata(my.recovery)



# Depth functions and descriptive analysis
# Exploratory FDA
# FPCA

# = =
id <- 1:261
dataf <- as.data.frame(df_wildfires$CA3989212363920080621)
ldata <- list(df = dataf,
              X.acres = unlist(list_acres, use.names=FALSE),
              X.lc_ent = unlist(list_landcover_entropy, use.names=FALSE))
f2 <- df_wildfires$CA3989212363920080621 ~ X.acres + X.lc_ent
basis1 <- create.fourier.basis(rangeval = c(1, 139), 
                               nbasis = 5)
basis.x1 <- list(X.acres = basis1)
basis2 <- create.bspline.basis(rangeval = c(1, 139), 
                               nbasis = 5)
basis.b1 <- list(X.acres = basis2)
res.lm2 <- fregre.lm(formula = f2, 
                     data = dataf, 
                     basis.x = basis.x1, 
                     basis.b = basis.b1)



ind <- 1:165
dataf <- as.data.frame(tecator$y[ind, ])
newdataf <- as.data.frame(tecator$y[-ind, ])
ldata <- list(df = dataf, X = X, X.d1 = X.d1, X.d2 = X.d2)
f2 <- Fat ~ Water + X.d2
basis.x1 <- list(X.d2 = basis1)
basis2 <- create.bspline.basis(rangeval = rangett, nbasis = 5)
basis.b1 <- list(X.d2 = basis2)
res.lm2 <- fregre.lm(f2, ldata, basis.x = basis.x1, basis.b = basis.b1)


absorp <- tecator$absorp.fdata

ind <- 1:165
tt <- absorp[["argvals"]]
y <- tecator$y$Fat[ind]
X <- absorp[ind, ]
X.d1 <- fdata.deriv(X, nbasis = 19, nderiv = 1) R> 
X.d2 <- fdata.deriv(X, nbasis = 19, nderiv = 2)

# = =   
  
class(my.recovery)
plot(my.recovery, axes=FALSE)
axis(1, at=seq(0, 139, 26), labels=c(seq(0, 7)))

plot(my.recovery, col='grey')
a1<-func.mean(my.recovery)
abline(h = mean(my.recovery$data[,1:125]), col='red')
abline(h = mean(my.recovery$data[,125:250]), col='green')
abline(v = my.recovery$argvals[125], col='black')
lines(a1, lwd=2)
points(a1$argvals, a1$data)

# How can we create a fd object (fda library) from the original elements?

?Data2fd
rownames(df_wildfires) <- years_argvals
my.recovery.fd <- Data2fd(argvals=years_argvals, 
                          y=as.matrix(df_wildfires),
                          covariates = #rbind(matrix(as.numeric(list_years)), 
                                             matrix(as.numeric(list_acres)),
                          fdname=list(xlab='Years since wildfire',
                                      ylab='NDVI',
                                      main='Landsat Surface Reflectance'))

plot(my.recovery.fd, ylab='NDVI')

### EXPLORE DERIVATIVES

plot(fdata.deriv(my.recovery.fd, 1), main='Derivative of recoveries', ylim=c(-.5, .5))
plot(fdata.deriv(my.recovery.fd,2), main='Second Derivative of recoveries', ylim=c(-5, 5))
plot(fdata.deriv(my.recovery.fd,3), main='Second Derivative of recoveries')

# Fitting functional regression -> Scalar-on-function regression SoFRe

library(refund)
data("DTI")
FA.cca <- DTI[complete.cases(DTI$cca),]
FA.cca$ID <- factor(FA.cca$ID)

pfr(data=my.recovery)

?pfr
sofr.fit <- pfr(pasat ~ lf(cca, k=30, argvals=1:93) + re(ID), data=FA.cca)
plot(sofr.fit, select=1, ylab=expression(paste(beta(t))), xlab="t (position along corpus callosum)")

# A function-on-scalar regression model:

fosr.fit <- pffr(cca~rbind(data.frame(list_years), data.frame(list_acres)), data=df_w)
plot(fosr.fit, pages=1, scale=0)

pffr(my.recovery.fd)

# FPCA

FPC.fit <- fpca.sc(FA.cca$cca, npc=10)
FPC.fit <- fpca.sc(my.recovery$data, npc=10)
cumsum(FPC.fit$evalues) / sum(FPC.fit$evalues) # First two values explain

pairs(FPC.fit$scores[,1:4]) # first and second component show that there is no non-linearities

df_lm <- df_covariates[df_covariates$Fire_ID %in% names(df_trends),]
df_lm$Y_scores <- FPC.fit$scores[,2]
lr <- lm(Y_scores ~ ., data=df_lm[, c('Year', 'Y_scores', 'StartMonth', 'StartDay', 'Acres',
                                      'Fire_Type', 'landcover', 'entropy_la')])
summary(lr)

pffr(Y_scores ~ c(s(Fire_Type)), data=df_lm, )

df_lm$response <- my.recovery$data

sofr <- pffr(response ~ Acres + Year + landcover + StartMonth + entropy_la, data=df_lm[df_lm$landcover != 'developed_medium_intensity',])
summary(sofr)
plot(sofr, pages=1)

df_lm$Year
df_lm$Acres <- scale(df_lm$Acres)


sofr_year <- pffr(response ~ c(Year) + Acres + StartMonth + entropy_la, yind=years_argvals, data=df_lm)
summary(sofr_year)
plot(sofr_year, pages=1)

# edf -> equivalent degrees of freedom, si es 1, es una constante

# s specifies non-parametric parameter
# c indicates the parameter to affect constantly
# re(ID) son los random effects





