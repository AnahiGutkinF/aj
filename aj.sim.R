####===========================================================================
# [aj] atkinson and juola recognition memory model
# Simulation Script
# Author: Anahi Gutkin
# Date: 17/09/2024
# Last Update: v5 (16/12/2024)
#####===========================================================================

#-------------------------------------------------------------------------------
# 1. Load packages and functions
#-------------------------------------------------------------------------------

# Go to the menu Session > Set Working Directory > To Source GitHub
# install.packages("methods")
# install.packages("numDeriv", repos = "http://cran.us.r-project.org")
# install.packages("matrixStats", repos = "http://cran.us.r-project.org")
# install.packages("statmod", repos = "http://cran.us.r-project.org")
# install.packages("devtools")
# install.packages("ggplot2")
# install.packages("patchwork")
# install.packages("tidyr")
library(tidyr)
library(MPTinR)
library(numDeriv)
library(doParallel)
library(parallel)
library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(Rmisc)
library(patchwork)
library(gridExtra)
library(data.table)

source("aj.fun.R")

#-------------------------------------------------------------------------------
# 2. Specify study conditions (Study Design)
#-------------------------------------------------------------------------------

#-----
# 2.1. Simulation conditions
#-----

r <- 1:2000 # Replications

n <- 100 # Sample size per tree.
g.m <- c("2ht", "sdt", "aj", "aj.3clk", "aj.3gk",
         "2ht0", "sdt0") # Generation model
par <- list( # Parameter
            c("dn" = 0.3, "do" = 0.5,
              "g3" = 0.4,
              "Hdnf" = 0.6, "Hdns" = 0.6,
              "Hdof" = 0.7, "Hdos" = 0.5, 
              "Ldn" = 0.7, "Ldo" = 0.9,
              "s_g1"= 0.7, "s_g2"= 0.8,
              "s_Hgnf"= 0.8, "s_Hgns"= 0.6, "s_Hgof"= 0.8, "s_Hgos"=0.6,
              "s_Lgn"= 0.6,  "s_Lgo"= 0.5),
            

            c("cr11" = -0.7, "cr12" = -0.5, "cr13" = -0.3,
              "cr2" = 1, "cr3" = 0.5,
              "d" = 1.5,
              "Lh" = 0.6, "Lr" = 0.5,
              "s_Lf"=0.8, "s_Lm"=0.7 ,
              "ss" = 1.5),
            
            c("c" = 1.5,
              "cl1" = -0.7, "cl2" = -0.4, "cl3" = -0.2,
              "d" = 1.5, "g3" = 0.4,
              "Hsf" = 0.9, "Hss" =  0.7,
              "Lf" = 0.8,    
              "r" = 0.6,   
              "s_g1" = 0.8, "s_g2" = 0.9, "s_Hff" = 0.9,
              "s_Hfs" = 0.8, "s_Hgf" = 0.6,"s_Hgs" = 0.5,
              "s_Lgn" = 0.5,"s_Lgo" = 0.6,"s_Ls" = 0.7,  
              "ss" = 0.7),
            
            c("c" = 1.5,
              "cl1" = -0.7, "cl2" = -0.4, "cl3" = -0.2,
              "d" = 1.5, "g" = 0.4,
              "Hsf" = 0.9, "Hss" =  0.7,
              "Lf" = 0.8,
              "r" = 0.6,
              "s_Hff" = 0.9,"s_Hfs" = 0.8, "s_Hgf" = 0.6,"s_Hgs" = 0.5,
              "s_Lgn" = 0.5,"s_Lgo" = 0.6,"s_Ls" = 0.8,
              "ss" = 1),
            
            c("c" = 1.5,
              "cl" = 1,
              "d" = 0.4, "g3" = 0.3,
              "Hsf" = 0.5, "Hss" = 0.5,
              "Lf" = 0.6,
              "r" = 0.4,
              "s_g1" = 0.8, "s_g2"=0.7,
              "s_Hff" = 0.9, "s_Hfs" = 0.8, "s_Hgf" = 0.4, "s_Hgs" = 0.6,
              "s_Lgn" = 0.5, "s_Lgo" = 0.6, "s_Ls" = 0.7,
              "ss" = 0.5),


            c("g1" = 0.5, "g2" = 0.4, "g3" = 0.3,
              "Hn" = 0.4, "Ho" = 0.5,
              "Ln" = 0.6, "Lo" = 0.7),
            
            c("cr11" = -0.7, "cr12" = -0.4, "cr13" = -0.2,
              "cr2" = 0.4, "cr3" = 0.5,
              "Hn" = 0.4, "Ho" = 0.5)
            )

# Create tibble
g.conds <- expand_grid(tibble(g.m = g.m, n = n, par = par), r = r)

#-----
# 2.2. Fitting conditions
#-----

f.m <- c("2ht", "sdt", "aj", "aj.3clk", "aj.3gk",
         "2ht0", "sdt0")

lb <- list(setNames(rep(0,17),
                    names(par[[1]])) ,
           setNames(c(rep(-Inf,3), 0,0, 0.1, rep(0,4), 0.1),
                    names(par[[2]])),
           setNames(c(0.1, rep(-Inf,3), 0.1, rep(0,14), 0.1),
                    names(par[[3]])),
           setNames(c(0.1, rep(-8, 3), 0.1, rep(0,12), 0.1),
                    names(par[[4]])),
           setNames(c(0.1, -8, 0.1, rep(0,14), 0.1),
                    names(par[[5]])),
           setNames(rep(0,7),
                    names(par[[6]])),
           setNames(c(rep(-Inf,3),0,0, 0, 0),
                    names(par[[7]]))
           )

hb <- list(setNames(rep(1,17),
                    names(par[[1]])) ,
           setNames(c(rep(Inf, 6), rep(1, 4), Inf),
                    names(par[[2]])),
           setNames(c(rep(Inf,5), rep(1,14), Inf),
                    names(par[[3]])),
           setNames(c(rep(Inf, 5), rep(1, 12), Inf),
                    names(par[[4]])),
           setNames(c(rep(Inf, 3), rep(1, 14), Inf),
                    names(par[[5]])),
           setNames(rep(1,7),
                    names(par[[6]])),
           setNames(c(rep(Inf,5), 1, 1),
                    names(par[[7]]))
           )

f.conds <- expand_grid(g.conds, tibble(f.m=f.m, lb=lb, hb=hb))

#-------------------------------------------------------------------------------
# 3. Simulation study
#-------------------------------------------------------------------------------

# Initial time
ini_time <- Sys.time()

#-----
# 3.a. Prepare for parallel
#-----

n.cores <- parallel::detectCores(logical = FALSE)-1
cl <- parallel::makeCluster(n.cores, type = "SOCK") 
doParallel::registerDoParallel(cl, cores = n.cores)
#-----
# 3.b. Simulation loop 
#-----

pm <- foreach::foreach(i = 1:nrow(f.conds), # performance metrics (pm)
                       .packages = c("MPTinR", "dplyr", "tidyr"), 
                       .combine = bind_rows, 
                       .inorder = FALSE
) %dopar% {
  try({ 
    
    #-----
    # 3.1. Establish randomization seed
    #-----
    
    set.seed(1320 + i) # Unique and identifiable seed for each condition
    
    #-----
    # 3.2. Specify condition
    #-----
    r <- f.conds$r[i]
    n <- f.conds$n[i]
    par <- f.conds$par[[i]]
    lb <- f.conds$lb[[i]]
    hb <- f.conds$hb[[i]]
    
    g.m.file <- paste0(getwd(), "/model_files/m.",
                       f.conds$g.m[i], ".txt") # Data-generating model file
    
    f.m.file <- paste0(getwd(), "/model_files/m.",
                       f.conds$f.m[i], ".txt")
    
    #-----
    # 3.3. Generate data
    #-----
    
    data <- MPTinR::gen.data(parameter.values = par, samples = 1,
                             model.filename = g.m.file,
                             n.per.item.type = rep(n,6), check.model = T)

    # Save Data (Optional)
    
    #-----
    # 3.4. Fit models and calculate performance measures
    #-----
    
    fit <- MPTinR::fit.model(data, model.filename = f.m.file, ci = .95,
                             lower.bound = lb ,upper.bound = hb)
    
    estim <- setNames(fit$parameters$estimates,
                      row.names(fit$parameters))
    
    my_list <- list(gof = fit$goodness.of.fit,
                    info = as.data.frame(fit$information.criteria),
                    estim=as.data.frame(t(estim)) )
    
    my_tibble <-f.conds[i,]
    for (name_list in names(my_list)) {
      my_tibble[[1, name_list]] <- my_list[[name_list]]
    }
    
  my_tibble # Final output tibbles
  })
}

parallel::stopCluster(cl) # Shutting down parallel processors

# Compute time 
fin_time <- Sys.time(); 
print(elapse_time <- fin_time - ini_time)

#-------------------------------------------------------------------------------
# 4. Save performance Measure
#-------------------------------------------------------------------------------

save(pm, file = "pm.R")
load("pm.R")

filter_pm <- function(data, f_m_filter, g_m_filter) {
  data %>%
    filter(f.m %in% f_m_filter, g.m %in% g_m_filter)
}

pm_1 <- filter_pm(pm, c("sdt", "2ht", "aj"),
                  c("sdt", "2ht", "aj", "sdt0", "2ht0"))

#-------------------------------------------------------------------------------
# 4. Performance Measures
#-------------------------------------------------------------------------------
#-----
# 4.1. Calculate nml
#-----

devtools::install_github("davidkellen/NMLmulti")
library(NMLmulti)

Rcpp::Rcpp.package.skeleton("modelsdt", attributes = TRUE, 
                            cpp_files = "modelsdt.cpp")
Rcpp::Rcpp.package.skeleton("model2ht", attributes = TRUE, 
                            cpp_files = "model2ht.cpp")
Rcpp::Rcpp.package.skeleton("modelaj", attributes = TRUE, 
                            cpp_files = "modelaj.cpp")
devtools::install("modelsdt")
devtools::install("model2ht")
devtools::install("modelaj")

library("model2ht")
library("modelaj")
library("modelsdt")

nml_sdt <- run_nml(fun= MSDT, parl=11, ks=rep(8,6), Ns=c(rep(50,6)),
                   packages_multicore = "modelsdt",
                   fits = 2, batchsize=1000, burn=10000, precision=.141,
                   cores = 20)

nml_2ht <- run_nml(fun= M2HT, parl=17, ks=rep(8,6), Ns=c(rep(50,6)),
                   packages_multicore = "model2ht",
                   fits = 2, batchsize=1000, burn=10000, precision=.141,
                   cores = 20)

nml_aj <- run_nml(fun= M2HT, parl=20, ks=rep(8,6), Ns=c(rep(50,6)),
                  packages_multicore = "modelaj",
                  fits = 2, batchsize=1000, burn=10000, precision=.141,
                  cores = 20)


save(nml_sdt,file="nmlsdt.RData")
save(nml_2ht,file="nml2ht.RData")
save(nml_aj,file="nmlaj.RData")

load("nmlsdt.RData")
load("nml2ht.RData")
load("nmlaj.RData")

nml <- list(nml_sdt, nml_2ht, nml_aj)

save(nml, file = "nml.RData")
nml_sdt # -68.13776
nml_2ht # -61.9646
nml_aj  # -63.40471

#-----
# 4.1. Calculate Lnml
#-----
pm_1 <- pm_1 %>%
  mutate(
    Lnml= case_when(
      f.m == "2ht" ~ gof$Log.Likelihood + nml_2ht$penalty,
      f.m == "sdt" ~ gof$Log.Likelihood + nml_sdt$penalty,
      f.m == "aj"  ~ gof$Log.Likelihood + nml_aj$penalty,
      TRUE ~ NA_real_  
    ),
    Ls = (
      gof$G.Squared/2 + gof$Log.Likelihood
    ),
    G2nml = case_when(
      f.m == "2ht" ~ 2*(Ls-Lnml),
      f.m == "sdt" ~ 2*(Ls-Lnml),
      f.m == "aj"  ~ 2*(Ls-Lnml),
      TRUE ~ NA_real_
    ),
    AICnml = case_when(
      f.m == "2ht" ~ G2nml+ 2*17,
      f.m == "sdt" ~ G2nml + 2*11,
      f.m == "aj"  ~ G2nml + 2*20,
      TRUE ~ NA_real_  
    )
  ); print(pm_1)


#-------------------------------------------------------------------------------
# 5. Model Recovery
#-------------------------------------------------------------------------------
# install.packages("dplyr")
library(dplyr)
library(writexl)

get_min_criteria <- function(data, criteria) {
  result <- data %>%
    group_by(r, g.m) %>%
    filter(info[[criteria]] == max(info[[criteria]], na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(g.m, f.m) %>%
    summarize(!!paste0("min.", criteria) := n(), .groups = "drop")
  return(result)
}

sel.AIC <- pm_1 %>%
                group_by(r, g.m) %>%
                filter(info$AIC == min(info$AIC, na.rm = TRUE)) %>%
                ungroup() %>%
                group_by(g.m, f.m) %>%
                tally(name = "min.AIC")

sel.BIC <- pm_1 %>%
                group_by(r, g.m) %>%
                filter(info$BIC == min(info$BIC, na.rm = TRUE)) %>%
                ungroup() %>%
                group_by(g.m, f.m) %>%
                tally(name = "min.BIC")

sel.Lnml <- pm_1 %>%
                group_by(r, g.m) %>%
                filter(Lnml == max(Lnml, na.rm = TRUE)) %>%
                ungroup() %>%
                group_by(g.m, f.m) %>%
                tally(name = "max.Lnml")

sel.G2nml <- pm_1 %>%
                group_by(r, g.m) %>%
                filter(G2nml == min(G2nml, na.rm = TRUE)) %>%
                ungroup() %>%
                group_by(g.m, f.m) %>%
                tally(name = "min.G2nml")

sel.AICnml <- pm_1 %>%
                group_by(r, g.m) %>%
                filter(AICnml == min(AICnml, na.rm = TRUE)) %>%
                ungroup() %>%
                group_by(g.m, f.m) %>%
                tally(name = "min.AICnml")

sel.Lnml <- mutate(sel.Lnml, "max%Lnml"=max.Lnml/2000)
sel.G2nml <- mutate(sel.G2nml, "min%G2nml"=min.G2nml/2000)
sel.AICnml <- mutate(sel.AICnml, "min%AICnml"=min.AICnml/2000)
sel.AIC <- mutate(sel.AIC, "min%AIC"=min.AIC/2000)
sel.BIC <- mutate(sel.BIC, "min%BIC"=min.BIC/2000)

#-------------------------------------------------------------------------------
# 6. Goodness of Fit
#-------------------------------------------------------------------------------

(gof_1 <- pm_1 %>%
  dplyr::group_by(g.m, f.m, .drop = FALSE) %>% 
  dplyr::summarize(RRate.GoF = round(mean(gof$p.value < 0.05, na.rm = TRUE),2),
                   .groups = "drop"))

(gof_2 <- pm_2 %>%
  dplyr::group_by(g.m, f.m, .drop = FALSE) %>%  
  dplyr::summarize(RRate.GoF = round(mean(gof$p.value < 0.05, na.rm = TRUE),2),
                   .groups = "drop"))
#-------------------------------------------------------------------------------
# 7. Merge performance measures
#-------------------------------------------------------------------------------

select <- gof_1 %>%
  full_join(sel.Lnml, by = c("g.m", "f.m"))%>%
  full_join(sel.AICnml, by = c("g.m", "f.m")) %>%
  full_join(sel.AIC, by = c("g.m", "f.m")) %>%
  full_join(sel.BIC, by = c("g.m", "f.m"))

writexl::write_xlsx(select, path = "select.xlsx")
