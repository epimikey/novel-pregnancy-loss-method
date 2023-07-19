###********************************************************************
###
### R code for analysis in:
###
### "Leung M, Rowland ST et al. A novel approach for inferring effects
###  of prenatal exposures on pregnancy loss"
###
###********************************************************************

###*************
### N: Notes ###
###*************

#' @param study_week study week number
#' @param week_start date of the beginning of the study week
#' @param lbic the number of LBIC for that week
#' @param avg_no2 average NO2 for that week
#' @param avg_temp average temperature for that week
#' @param doy day of the year of the beginning of the study week
#' @param year year of the study week
#' @param month month of the study week
#' @param season season of the study week
#' @param harm within-year harmonic of the beginning of the study week
#' @param ndlag lag of weekly NO2 (e.g., ndlag01 is week 1 lag of NO2, ndlag02 is week 2 lag of NO2 etc.)
#' @param templag lag of weekly temperature (e.g., templag01 is week 1 lag of temperature, templag02 is week 2 lag of temperature)

###*************************
### 1: Required Packages ###
###*************************

# install required packages
requiredPackages = c('tidyverse','dlnm','splines','boot','ggplot2')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

####*********************************************
#### 2: Setup the data for DLM g-computation ####
####*********************************************

# load dataset
df <- readRDS("lbic_mock_df.rds")

# create exposure histories and cross basis
nd.hist.obs <- df %>% select(ndlag01:ndlag40)
temp.hist.obs <- df %>% select(templag01:templag40)

# test AIC for model parameterization
dfAIC <- data.frame(df1 = numeric(), df2 = numeric(), aic = numeric())
for (i in 4:10) {
  cb.nd.obs <- crossbasis(nd.hist.obs, argvar = list("lin"), arglag = list(df=i))
  for (j in 4:10) {
    cb.temp.obs <- crossbasis(temp.hist.obs, argvar = list(df=3), arglag = list(df=j))
    modelAIC <- glm(lbic ~ cb.nd.obs + cb.temp.obs + as.factor(year)*harm, family=poisson(), data=df)
    dfAIC <- bind_rows(dfAIC, data.frame(df1 = i, df2 = j, aic = AIC(modelAIC)))
  }
}
dfAIC <- dfAIC %>% filter(aic == min(aic)) # 5 df for NO2, 4 df for temperature

# define crossbasis and merge with original df
cb.nd.obs <- crossbasis(nd.hist.obs, argvar = list("lin"), arglag = list(df=dfAIC[1,1])) %>% data.frame()
cb.temp.obs <- crossbasis(temp.hist.obs, argvar = list(df=3), arglag = list(df=dfAIC[1,2])) %>% data.frame()
names(cb.nd.obs) <- paste0('cb.', "nd", '.', names(cb.nd.obs))
names(cb.temp.obs) <- paste0('cb.', "temp", '.', names(cb.temp.obs))
df <- df %>% bind_cols(cb.nd.obs, cb.temp.obs)

####***************************************************************
#### 3: Write functions to generate counterfactual predictions ####
####***************************************************************

# create counterfactual exposure histories and crossbasis
make_cfNO2 <- function(cfNO2) {
  
  # make a dataframe where each column one gestational week has NO2 of 
  # the set value and NA otherwise
  cfNO2.df <- diag(cfNO2, nrow = 40, ncol = 40)
  cfNO2.df <- as.data.frame(cfNO2.df)
  cfNO2.df[cfNO2.df==0] <- NA
  
  # add a column for when all weeks are exposed
  cfNO2.df[,41] <- rep(cfNO2, 40)
  
  # transpose the dataframe so we can then split the columns into a list
  cfNO2.df <- t(cfNO2.df) %>% 
    as.data.frame() %>% 
    dplyr::mutate(expProfile = row_number())
  
  # remove row names
  rownames(cfNO2.df) <- NULL
  
  # return the list of exposure profiles
  return(cfNO2.df)
}

cf_cb <- function(expdf, covariate, cfval, w) {
  
  # create cf exposure history
  cfNO2 <- make_cfNO2(cfval)
  
  # identify column to replace
  if(w < 41 & !is.na(cfNO2[w,w])){
    # name of column we will replace
    varName.no2 <- paste0("ndlag", stringr::str_pad(w, 2, "left", "0" ))
    # now replace that column with the counterfactual exposure value
    expdf[ ,varName.no2] <- cfNO2[w,w]
  }
  if(w ==41 & !is.na(cfNO2[w,w])){
    expdf[,] <- cfval
  }
  
  cb <- crossbasis(expdf, argvar = list("lin"), arglag = list(df=dfAIC[1,1])) %>% data.frame()
  names(cb) <- paste0('cb.', "nd", '.', names(cb))
  cb <- cbind(cb,covariate) %>% mutate(cfexp = cfval, week = w)
  cb
}

####****************************************************
#### 4: Count Difference - Pregnancy LOss (CD_{PL}) ####
####****************************************************

lbic.std.boot <- function(data, indices) {
  
  # display progress with each run of the function
  d <- data[indices,]
  p$tick()$print()  # update progress bar
  
  # exposure matrix
  nd.hist <- d %>% select(ndlag01:ndlag40)
  
  # isolate covariates
  covariate <- d %>% select(cb.temp.v1.l1:cb.temp.v3.l4, year, harm)
  
  # fit Q model (i.e., using observed data)
  Qeq <- glm(lbic ~ cb.nd.v1.l1 + cb.nd.v1.l2 + cb.nd.v1.l3 + cb.nd.v1.l4 + cb.nd.v1.l5 + 
               cb.temp.v1.l1 + cb.temp.v1.l2 + cb.temp.v1.l3 + cb.temp.v1.l4 + 
               cb.temp.v2.l1 + cb.temp.v2.l2 + cb.temp.v2.l3 + cb.temp.v2.l4 + 
               cb.temp.v3.l1 + cb.temp.v3.l2 + cb.temp.v3.l3 + cb.temp.v3.l4 + 
               as.factor(year)*harm,
             family = poisson(link="log"), 
             data = d)
  
  # create empty df for results
  all_cb_ref <- d %>% select(cb.nd.v1.l1:cb.temp.v3.l4, year, harm) %>% filter(row_number() < 1)
  
  # create counterfactual exposure and exposure histories
  for (i in 1:2) {
    # loop through exposure histories: 0, 10, 20, 30, 40 ppb
    cfval <- i*10
    # loop through gestational weeks
    for (j in 1:41) {
      all_cb_ref <- bind_rows(all_cb_ref, cf_cb(nd.hist, covariate, cfval, j))
    }
  }
  # predict and take the average for each bootstrap iteration
  all_cb_ref$lbic_pr <- predict(Qeq, all_cb_ref, type="response", se.fit=F)
  pr <- all_cb_ref %>% 
    group_by(cfexp, week) %>%
    summarise(mean_lbic_pr = mean(lbic_pr)) %>%
    mutate(cfexp = paste0("cfexp",cfexp)) %>%
    spread(cfexp, mean_lbic_pr) 
  
  # calculate contrasts
  pr <- pr %>% mutate(est = cfexp20 - cfexp10)
  pr_cum <- pr %>% 
    filter(week >= 1 & week <= 40) 
  return(c(-pr$est))
}

tot_rep <- 500 # number of bootstrap samples
p <- progress_estimated(tot_rep+1)
lbic.std.results <- boot(data = df, statistic = lbic.std.boot, R = tot_rep)

# extract the bootstrap results
lbic.std.boot.est <- data.frame(cbind(original=lbic.std.results$t0[1:41], t(lbic.std.results$t[,1:41])))

# calculate the mean and percentile-based 95% confidence intervals
boot.results.est <- data.frame(cbind(week=seq(1,41), est=lbic.std.boot.est$original),
                               lci=(apply((lbic.std.boot.est)[,-1], 1, quantile, probs=0.025)),
                               uci=(apply((lbic.std.boot.est)[,-1], 1, quantile, probs=0.975)))

# from boot.results.est, week 41 gives you the cumulative estimate over the lag period


####*******************************************************************
#### 5. Risk Ratio - Live birth identified conceptions (RR_{LBIC}) ####
####*******************************************************************

nd.hist.obs <- df %>% select(ndlag01:ndlag40)
temp.hist.obs <- df %>% select(templag01:templag40)

cb.nd.obs <- crossbasis(nd.hist.obs, argvar = list("lin"), arglag = list(df = dfAIC[1,1]))
cb.temp.obs <- crossbasis(temp.hist.obs, argvar = list(df = 3), arglag = list(df = dfAIC[1,2]))

RRmodel <- glm(lbic ~ cb.nd.obs + cb.temp.obs + as.factor(year) + harm,
               family = quasipoisson(link = "log"),
               data = df)


pred.RRLBIC <- crosspred(cb.nd.obs, RRmodel, at=20, cen=10)

results.RRLBIC <- data.frame(
  week = 1:40,
  est = pred.RRLBIC$matRRfit[1,],
  lci = pred.RRLBIC$matRRlow[1,],
  uci = pred.RRLBIC$matRRhigh[1,],
  type = 2
)

####***********************************************
#### 6. Plotting CD_{PL} and RR_{LBIC} Results ####
####***********************************************

# combine the CD_{PL} and RR_{LBIC} results
results <- boot.results.est %>% filter(week < 41) %>% mutate(type = 1) %>%
  bind_rows(results.RRLBIC) %>%
  mutate(type = factor(type,
                       levels=c(1,2),
                       labels=c(expression(paste("CD"[PL])),expression(paste("RR"[LBIC])))))

# create the different y-intercepts (0 for CD_{PL} and 1 for RR_{LBIC})
yint <- data.frame(type = 1:2, int = 0:1) %>% 
  mutate(type = factor(type,
                       levels=c(1,2),
                       labels=c(expression(paste("CD"[PL])),expression(paste("RR"[LBIC])))))

# plot the results
results %>%
  ggplot() +
  geom_ribbon(aes(x=week, ymin=uci, ymax=lci), fill="grey", alpha=0.7) +
  geom_line(aes(x=week, y=est), color="black") +
  geom_hline(data=yint, aes(yintercept=int), linetype=2, color="black") +
  labs(x = "Gestational week", y = expression(paste("Estimate per 10 ppb NO"[2]))) +
  facet_grid(type~., scales="free_y", labeller = label_parsed) +
  scale_x_continuous(breaks=seq(-15, 40, 5)) + 
  theme(panel.spacing.y=unit(1.5,"lines"),
        panel.spacing.x=unit(1.25,"lines"),
        panel.border=element_rect(colour="black",size=0.2,fill=NA),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.y=element_text(size=16,angle=0, hjust=0))


