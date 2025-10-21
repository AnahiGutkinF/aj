#####===========================================================================
##### [aj.mpt] A-J multinomial process tree model
##### Fit Data Script
##### Autor: Anahi Gutkin
##### Date: 27/03/2024
##### Version: v5 (14/10/2025)
#####===========================================================================

#-------------------------------------------------------------------------------
# 1. Set conditions, Load packages, functions, data
#-------------------------------------------------------------------------------

# set working directory

packages <- c("MPTinR", "openxlsx", "snow", "ggplot2", "tidyr", "patchwork",
              "data.table", "parallel", "writexl", "caret", "snowfall")
lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

# source("aj.model.R")
source("aj.fun.R")

# -------------------------------------------------------------------------------
# 2. Conditions & Data
#-------------------------------------------------------------------------------

objetos <- ls()
cond <- c("aj.3clk","aj.3gk", "sdt", "2ht")
f <- paste0("fit.", cond) 
m <- paste0("m.", cond) 
conditions <- cbind("fit" = f, "model"=m); rm(cond, f, m, objetos)
load("d.data.cl.rt.RData")
n_sub <- nrow(d.data.cl.rt)# Number of subjects
n.cores <- parallel::detectCores()-1


#-------------------------------------------------------------------------------
# 3. Main Models
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 3.1 Check
#-------------------------------------------------------------------------------

m.2ht <- paste0(getwd(), "/model_files", "/m.2ht.txt")
m.sdt <- paste0(getwd(), "/model_files", "/m.sdt.txt")
m.aj <- paste0(getwd(), "/model_files", "/m.aj.txt")
m.2ht.5k <- paste0(getwd(), "/model_files", "/m.2ht.5k.txt")
m.sdt.5k <- paste0(getwd(), "/model_files", "/m.sdt.5k.txt")
m.aj.5k <- paste0(getwd(), "/model_files", "/m.aj.5k.txt")

check.mpt(m.2ht)
check.mpt(m.sdt)
check.mpt(m.aj)

check.mpt(m.2ht.5k)
check.mpt(m.sdt.5k)
check.mpt(m.aj.5k)

#-------------------------------------------------------------------------------
# 3.2 Fit
#-------------------------------------------------------------------------------

fit.2ht<-  fit.model(model.filename = m.2ht, data = d.data.cl.rt, 
                     multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                     n.optim =  20,
                     ci = 95, 
                     use.gradient = F);
save(fit.2ht, file = "fit.2ht.RData") 

fit.sdt<-  fit.model(model.filename = m.sdt, data = d.data.cl.rt, 
                     multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                     n.optim =  20,
                     ci = 95, 
                     lower.bound= c(-Inf, 0, 0, 0.1, rep(0,6), 0.1),
                     upper.bound= c(rep(Inf, 4), rep(1, 2), rep(Inf,2), rep(1,2), Inf),
                     use.gradient = F);save(fit.sdt, file = "fit.sdt.RData")

fit.aj <-   fit.model(model.filename = m.aj, data = d.data.cl.rt, 
                      multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                      n.optim =  20,
                      ci = 95, 
                      output = "full",
                      lower.bound= c(0.1, rep(-Inf,1), 0.1, rep(0,16), 0.1),
                      upper.bound= c(rep(Inf,3), rep(1,5), Inf, Inf, rep(1,9), Inf)
                      );save(fit.aj, file = "fit.aj.RData")

# check.mpt(m.aj, restrictions.filename = list("s_g1=1", "s_g2=1"))
fit.aj.3clk <-  fit.model(model.filename = m.aj, data = d.data.cl.rt, 
                          multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          output = "full",
                          lower.bound= c(0.1, rep(-8, 1), 0.1, rep(0,14), 0.1),
                          upper.bound= c(rep(Inf,3), rep(1, 5), rep(Inf,2), rep(1,7), 1),
                          restrictions.filename = list("s_g1=1", "s_g2=1")
                          );save(fit.aj.3clk, file = "fit.aj.3clk.RData")

# check.mpt(m.aj, restrictions.filename = list("s_cl2=1", "s_cl3=1"))
fit.aj.3gk <-  fit.model(model.filename = m.aj.3gk, data = d.data.cl.rt, 
                         multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          output = "full",
                          lower.bound= c(0.1, -8, 0.1, rep(0,14), 0.1),
                          upper.bound= c(rep(Inf, 3), rep(1, 14), Inf),
                          restrictions.filename = list("s_cl2=1", "s_cl3=1"),
                          use.gradient = T);save(fit.aj.3gk, file = "fit.aj.3gk.RData")

load_files <- c("fit.2ht.RData", "fit.sdt.RData", "fit.aj.RData", 
                "fit.aj.3clk.RData", "fit.aj.3gk.RData")

for (file in load_files) {
  load(file = file)
}

#-------------------------------------------------------------------------------
# 3 Model comparisons
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 3.1 Non nested model comparisons: Main models
#-------------------------------------------------------------------------------

sel_list1 <- select.mpt(list(fit.2ht, fit.sdt, fit.aj))
sel_list2 <- select.mpt(list(fit.2ht, fit.sdt, fit.aj.3clk, fit.aj.3gk))

if (ncol(sel_list1) == ncol(sel_list2)) {
  sel <- rbind(sel_list1, sel_list2)
  
  sel <- transform(sel, # Table 1
                   "p%smaller.05" = (p.smaller.05 / 47) * 100, 
                   "AIC%best" = (AIC.best / 47) * 100)[, c("model",
                                                           "n.parameters",
                                                           "p.smaller.05",
                                                           "p%smaller.05",
                                                           "AIC.best",
                                                           "AIC%best")]
} else {
  stop("Error: list with different number of columns.")
}



format_p <- function(p_values) {
  sapply(p_values, function(p) {
    if (p < 0.001) {
      "<.001"
    } else {
      sub("^0", "", format(round(p, 3), nsmall = 3))
    }
  })
}

format_x2_aic <- function(values) {
  format(round(values, 2), nsmall = 2)
}


tabla1 <- data.frame(   # Table B.2
  Id = 1:47,
  `X2(df=51)` = format_x2_aic(fit.2ht$goodness.of.fit$individual$G.Squared),
  p_2HT = format_p(fit.2ht$goodness.of.fit$individual$p.value),
  AIC_2HT = format_x2_aic(fit.2ht$information.criteria$individual$AIC),
  
  `X2(df=57)` = format_x2_aic(fit.sdt$goodness.of.fit$individual$G.Squared),
  p_SDT = format_p(fit.sdt$goodness.of.fit$individual$p.value),
  AIC_SDT = format_x2_aic(fit.sdt$information.criteria$individual$AIC),
  
  `X2(df=50)_Ajcj` = format_x2_aic(fit.aj.3clk$goodness.of.fit$individual$G.Squared),
  p_Ajcj = format_p(fit.aj.3clk$goodness.of.fit$individual$p.value),
  AIC_Ajcj = format_x2_aic(fit.aj.3clk$information.criteria$individual$AIC),
  
  `X2(df=50)_AJgj` = format_x2_aic(fit.aj.3gk$goodness.of.fit$individual$G.Squared),
  p_AJgj = format_p(fit.aj.3gk$goodness.of.fit$individual$p.value),
  AIC_AJgj = format_x2_aic(fit.aj.3gk$information.criteria$individual$AIC)
)

# writexl::write_xlsx(tabla1, path = "tableb2.xlsx")
#-------------------------------------------------------------------------------
# 3.2 Nested model comparisons
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 3.2.1 Testing Predictions
#-------------------------------------------------------------------------------

check.mpt(m.aj, restrictions.filename = list("s_g1=1", "s_g2=1",
                                             "s_cl2=0", "s_cl3=0"))
fit.aj.3clk.r1 <- fit.model(model.filename = m.aj, data = d.data.cl.rt,
                            multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                            n.optim =  20,
                            ci = 95,
                            output = "full",
                            lower.bound= c(0, rep(-8, 1), 0.1, rep(0,12), 0.1),
                            upper.bound= c(rep(Inf,3), rep(1, 12), Inf),
                            restrictions.filename = list("s_g1=1", "s_g2=1",
                                                         "s_cl2=0", "s_cl3=0"))

# save(fit.aj.3clk.r1, file = "fit.aj.3clk.r1.RData")
load("fit.aj.3clk.r1.RData")

nested.cj <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.r1)
nested.gj <- test.nested(model.gen = fit.aj.3gk , model.rest = fit.aj.3clk.r1)

sum(nested.cj$p.val<0.05)/47
sum(nested.gj$p.val<0.05)/47

format_p <- function(p_values) {
  sapply(p_values, function(p) {
    if (p < 0.001) {
      "<.001"
    } else if (p == 1) {
      ".999"  
    } else {
      sub("^0", "", format(round(p, 3), nsmall = 3))
    }
  })
}

table_nested1 <- data.frame( # Table 2.
  Id = 1:47,
  
  `G2 (df=3)_AJcj` = nested.cj$G2,
  `p_AJcj` = format_p(nested.cj$p.val),
  
  `G2 (df=3)_AJgj` = nested.gj$G2,
  `p_AJgj` = format_p(nested.gj$p.val)
)


#-------------------------------------------------------------------------------
# 3.2.3 Testing A-J assumptions about RTs (L parameters)
#-------------------------------------------------------------------------------

check.mpt(m.aj, restrictions.filename = list("s_Lgn=1", "s_Lgo=1", "s_Ls=1", 
                                             "s_g1=1", "s_g2=1"))
fit.aj.3clk.rL <-  fit.model(model.filename = m.aj,
                            data = d.data.cl.rt, 
                            ci = 95, 
                            n.optim =  20,
                            lower.bound= c(0, -Inf, 0.1, rep(0,11), 0.1),
                            upper.bound= c(rep(Inf,3), rep(1,5),
                                           rep(Inf, 2), rep(1,4), Inf),
                            restrictions.filename = list("s_Lgn=1", "s_Lgo=1",
                                                          "s_Ls=1",
                                                          "s_g1=1", "s_g2=1"))
# save(fit.aj.3clk.rL, file = "fit.aj.3clk.rL.RData")
load("fit.aj.3clk.rL.RData")


check.mpt(m.aj, restrictions.filename = list("s_Lgn=1", "s_Lgo=1", "s_Ls=1",
                                             "s_cl2=0", "s_cl3=0"))
fit.aj.3g.rL <-  fit.model(model.filename = m.aj,
                           data = d.data.cl.rt, 
                           ci = 95, 
                           n.optim =  20,
                           lower.bound= c(0, -Inf, 0.1, rep(0,11), 0.1),
                           upper.bound= c(rep(Inf,3), rep(1,11), rep(Inf,1)),
                           restrictions.filename = list("s_Lgn=1", "s_Lgo=1",
                                                        "s_Ls=1",
                                                        "s_cl2=0", "s_cl3=0"))
# save(fit.aj.3g.rL, file = "fit.aj.3g.rL.RData")
load("fit.aj.3g.rL.RData")

nested.cj.L <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.rL)
nested.gj.L <- test.nested(model.gen = fit.aj.3gk , model.rest = fit.aj.3g.rL)

sum(nested.cj.L$p.val<.05)/47
sum(nested.gj.L$p.val<.05)/47

table_nested2 <- data.frame( # Table 3.
  Id = 1:47,
  
  `G2 (df=3)_AJcj` = nested.cj.L$G2,
  `p_AJcj` = format_p(nested.cj.L$p.val),
  
  `G2 (df=3)_AJgj` = nested.gj.L$G2,
  `p_AJgj` = format_p(nested.gj.L$p.val)
)

# write_xlsx(table_nested2, "table_nested2.xlsx")

#-------------------------------------------------------------------------------
# 6.2.4 Testing A-J assumptions about CLs (H parameters)
#-------------------------------------------------------------------------------

check.mpt(m.aj, restrictions.filename = list("s_Hff=1", "s_Hfs=1",
                                             "s_Hgf=1","s_Hgs=1",
                                             "s_g1=1",  "s_g2=1"))

fit.aj.3clk.rH <-  fit.model(model.filename = m.aj,
                             data = d.data.cl.rt, 
                             n.optim =  20,
                             ci = 95, 
                             lower.bound= c(0.1, -Inf, 0.1, rep(0,10),0.1),
                             upper.bound= c(rep(Inf,3), rep(1,5),
                                            rep(Inf,2), rep(1,3), Inf), 
                             restrictions.filename = list("s_Hff=1", "s_Hfs=1",
                                                          "s_Hgf=1","s_Hgs=1",
                                                          "s_g1=1",  "s_g2=1"));
# save(fit.aj.3clk.rH, file = "fit.aj.3clk.rH.RData");
load("fit.aj.3clk.rH.RData")

check.mpt(m.aj, restrictions.filename = list("s_Hff=1", "s_Hfs=1",
                                             "s_Hgf=1","s_Hgs=1",
                                             "s_cl2=0",  "s_cl3=0"))

fit.aj.3gk.rH <-  fit.model(model.filename = m.aj,
                            data = d.data.cl.rt,
                            ci = 95, 
                            n.optim =  20,
                            lower.bound= c(0.1, -Inf, 0.1, rep(0,10), 0.1),
                            upper.bound= c(rep(Inf,3), rep(1,10), Inf),
                            restrictions.filename = list("s_Hff=1", "s_Hfs=1",
                                                         "s_Hgf=1","s_Hgs=1",
                                                         "s_cl2=0",  "s_cl3=0"));
# save(fit.aj.3gk.rH, file = "fit.aj.3gk.rH.RData");
load("fit.aj.3gk.rH.RData")

nested.cj.H <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.rH)
nested.gj.H <- test.nested(model.gen = fit.aj.3gk, model.rest = fit.aj.3gk.rH )


sum(nested.cj.H$p.val<0.05)/47
sum(nested.gj.H$p.val<0.05)/47
    
    table_nested3 <- data.frame( # Table 4.
      Id = 1:47,
      
      `G2 (df=4)_AJcj` = nested.cj.H$G2,
      `p_AJcj` = format_p(nested.cj.H$p.val),
      
      `G2 (df=4)_AJgj` = nested.gj.H$G2,
      `p_AJgj` = format_p(nested.gj.H$p.val)
    )
    
    # write_xlsx(table_nested3, "table_nested3.xlsx") # Optional

#-------------------------------------------------------------------------------
# 6.3. p-values plots
#-------------------------------------------------------------------------------
    
  plot_aj_logp_histogram <- function(gj_obj, cj_obj, label = "AJ") {
    # Build group names
    name_g <- paste0(label, "g")
    name_c <- paste0(label, "c")
    
    # Ensure flat numeric vectors
    gj_vals <- as.numeric(unlist(gj_obj$p.val))
    cj_vals <- as.numeric(unlist(cj_obj$p.val))
    
    # Replace zeros to avoid log(0)
    gj_vals[gj_vals == 0] <- 0.001
    cj_vals[cj_vals == 0] <- 0.001
    
    # Create individual data frames
    gj_df <- data.frame(
      log_pval = log(gj_vals),
      test = name_g
    )
    cj_df <- data.frame(
      log_pval = log(cj_vals),
      test = name_c
    )
    
    # Combine both into one data frame
    plot_df <- dplyr::bind_rows(gj_df, cj_df)
    
    # Threshold in log scale
    threshold_log <- log(0.05)
    
    # Create faceted histogram with vertical line and label at y = 47
    ggplot2::ggplot(plot_df, ggplot2::aes(x = log_pval, fill = test)) +
      ggplot2::geom_histogram(alpha = 0.9, bins = 30, color = "black") +
      ggplot2::geom_vline(xintercept = threshold_log, linetype = "dashed", color = "black", linewidth = 0.6) +
      ggplot2::annotate("text",
                        x = threshold_log,
                        y = 40,
                        label = "α = 0.05",
                        vjust = -0.2, hjust = -0.1,
                        size = 3, color = "black") +
      ggplot2::facet_wrap(~ test, ncol = 1) +
      ggplot2::labs(
        x = "log(p-value)",
        y = "Frequency",
        fill = "Model"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 47)) +
      ggplot2::scale_fill_manual(values = setNames(c("gray30", "gray60"), c(name_c, name_g))) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        axis.title = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 9),
        legend.position = "none"  # optional: hide legend since facets label groups
      )
  }
  
  
  
  p_AJ_hist    <- plot_aj_logp_histogram(nested.gj, nested.cj, label = "AJ")
  p_AJ_hist_H  <- plot_aj_logp_histogram(nested.gj.H, nested.cj.H, label = "AJ")
  p_AJ_hist_L  <- plot_aj_logp_histogram(nested.gj.L, nested.cj.L, label = "AJ")
  
  
  ggsave("aj_hist_hypothesis_testing.png", plot = p_AJ_hist,
         width = 8, height = 4, dpi = 300, units = "in")
  
  ggsave("aj__hist_hypothesis_testing_H.png", plot = p_AJ_hist_H,
         width = 8, height = 4, dpi = 300, units = "in")
  
  ggsave("aj__hist_hypothesis_testing_L.png", plot = p_AJ_hist_L,
         width = 8, height = 4, dpi = 300, units = "in")
  
#-------------------------------------------------------------------------------
# 4. Import fits 
#-------------------------------------------------------------------------------

archivos_fit <- list.files(getwd(), pattern = "^fit\\.")
for (archivo in archivos_fit) {
  ruta_archivo <- file.path(getwd(), archivo)
  load(ruta_archivo, .GlobalEnv)  # Carga el objeto en el entorno global
  cat("Imported data:", archivo, "\n")
  }

load("~/R.Projects/aj/nmlaj.RData")
load("~/R.Projects/aj/nml2ht.RData")
load("~/R.Projects/aj/nmlsdt.RData")

#-------------------------------------------------------------------------------
# 7. Grouping par estimations & Mean par
#-------------------------------------------------------------------------------

par.3clk <- do.call(rbind, lapply(1:n_sub, function(i) fit.aj.3clk$parameters$individual[,"estimates",i]))
par.gk <- do.call(rbind, lapply(1:n_sub, function(i) fit.aj.3gk$parameters$individual[,"estimates",i]))
colnames(par.3clk) <- names(fit.aj.3clk$parameters$individual[,"estimates",1])
colnames(par.gk) <- names(fit.aj.3gk$parameters$individual[,"estimates",1])
mean_3clk <- colMeans(par.3clk[, ])
mean_gk <- colMeans(par.gk[, ])
# Mean parameters
print(mean_3clk)
print(mean_gk)

#-------------------------------------------------------------------------------
# 9. Restructure data with 5 target conditions 
#-------------------------------------------------------------------------------

# Define categories and confidence levels
cat <- c("hit", "miss", "fa", "cr")
cl <- c("high", "medium", "low")
k <- paste0("k", rep(1:3, each = 24))
l.x <- paste0(c("hit", "miss", "fa", "cr"), "-k", rep(1:5, each = 4))
l.cnf <- paste0(rep(l.x, each = 3), c("-high", "-medium", "-low"))
l.cnf2 <- paste0(rep(l.x, each = 2), c("-high", "-low"))
cat.order <- paste0(rep(cat, each = 3 * 2), "-", k, "-", rep(cl, each = 2), c(1, 0))

# Read and preprocess data
d.juola <- read.delim("data_Juola.txt")
d.juola$Stimulus <- change_factor_levels(as.factor(d.juola$Stimulus), c("new_item", "learned_item"))
d.juola$Type_of_trial <- change_factor_levels(as.factor(d.juola$Type_of_trial), c("hit", "miss", "fa", "cr"))
d.juola$Confidence_level <- change_factor_levels(as.factor(d.juola$Confidence_level), c("high", "medium", "low"))
d.juola$Target_frequency2 <- d.juola$Target_frequency

# Create a copy of the data with modified confidence levels
d.juola2 <- d.juola
d.juola2$Confidence_level[d.juola2$Confidence_level == "medium"] <- "low"

# With 3 confidence levels
d.juola <- data.frame(
  tree = d.juola$Stimulus,
  x = factor(paste0(d.juola$Type_of_trial, "-k", d.juola$Target_frequency2), levels = l.x),
  x_cnf = factor(paste0(d.juola$Type_of_trial, "-k", d.juola$Target_frequency2, "-", d.juola$Confidence_level), levels = l.cnf),
  y = d.juola$RT,
  id = d.juola$Id
)

# Data for 5 %targets conditions & 3 confidence levels
data.5j.3CL <- recateg.2bin(data = d.juola, NSets = 47, criteria = "geom", groupby = "x_cnf")

cat.5j.3CL <- paste(
  rep(cat, each = 4),
  rep(c("k1", "k2", "k3", "k4", "k5"), each = 16),
  paste0(rep(cl, 2), rep(0:1, each = 3)),
  sep = "-"
)

# With 2 confidence levels
d.juola2 <- data.frame(
  tree = d.juola2$Stimulus,
  x = factor(paste0(d.juola2$Type_of_trial, "-k", d.juola2$Target_frequency2), levels = l.x),
  x_cnf = factor(paste0(d.juola2$Type_of_trial, "-k", d.juola2$Target_frequency2, "-", d.juola2$Confidence_level), levels = l.cnf2),
  y = d.juola2$RT,
  id = d.juola2$Id
)

# Data for 5 %targets conditions & 2 confidence levelsx
data.5j.2CL <- recateg.2bin(data = d.juola2, NSets = 47, criteria = "geom", groupby = "x_cnf")

cat.5j.2CL <- paste(
  rep(cat, each = 4),
  rep(c("k1", "k2", "k3", "k4", "k5"), each = 16),
  paste0(rep(c("high", "low"), 2), rep(0:1, each = 3)),
  sep = "-"
)

# Change order of data
data.5j.3CL <- data.5j.3CL[, cat.5j.3CL]
data.5j.2CL <- data.5j.2CL[, cat.5j.2CL]

# Save data
save(data.5j.3CL, file = "data.5j.3CL.RData")
save(data.5j.2CL, file = "data.5j.2CL.RData")



#-------------------------------------------------------------------------------
# 10. Model fitting for 5 %target conditions
#-------------------------------------------------------------------------------

check.mpt(m.aj.5k, restrictions.filename = list("s_g1=1", "s_g2=1",
                                                "s_g3=1", "s_g4=1"))
check.mpt(m.aj.5k, restrictions.filename = list("s_cl2=0", "s_cl3=0",
                                                "s_cl4=0", "s_cl5=0"))

n.cores <- parallel::detectCores()-1

fit.2ht.5k <- fit.model(model.filename = m.2ht.5k, data = data.5j.2CL, 
                         multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                         n.optim =  20,
                         ci = 95, 
                         use.gradient = F,
                         )
save(fit.2ht.5k, file = "fit.2ht.5k.RData") 

fit.aj.5k.3clk <- fit.model(model.filename = m.aj.5k, data = data.5j.2CL, 
                           multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                           n.optim =  20,
                           ci = 95,  
                           lower.bound = c(0.1, -8, 0.1, rep(0, 16), 0.1),
                           upper.bound = c(rep(Inf, 3), rep(1, 5), rep(1,4), rep(1,7), Inf),
                           use.gradient = F,
                           restrictions.filename = list("s_g1=1", "s_g2=1",
                                                        "s_g3=1", "s_g4=1"),
                           );save(fit.aj.5k.3clk, file = "fit.aj.5k.3clk.RData") 

fit.aj.5k.3gk <- fit.model(model.filename = m.aj.5k, data = data.5j.2CL, 
                       multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                       n.optim =  20,
                       ci = 95,  
                       lower.bound = c(0.1, -8, 0.1, rep(0, 16), 0.1),
                       upper.bound = c(rep(Inf, 3), rep(1, 16), Inf),
                       use.gradient = F,
                       restrictions.filename = list("s_cl2=0", "s_cl3=0",
                                                    "s_cl4=0", "s_cl5=0"),
                       );save(fit.aj.5k.3gk, file = "fit.aj.5k.3gk") 

fit.sdt.5k <- fit.model(model.filename = m.sdt.5k, data = data.5j.2CL, 
                            multicore = "n.optim", sfInit = TRUE, ci = 95, 
                            lower.bound = c(-8, rep(0,2), 0.1, rep(0, 8), 0.1),
                            upper.bound = c(rep(Inf, 4), rep(1, 2), rep(Inf,4), rep(1,2), Inf),
                            use.gradient = FALSE
                            );save(fit.sdt.5k, file = "fit.sdt.5k.RData")

load_files <- c("fit.2ht.5k.RData", "fit.sdt.5k.RData",
                "fit.aj.5k.3gk.RData", "fit.aj.5k.3clk.RData")

#Appendix B (Table B.1)
par <- setNames(fit.2ht.5k$parameters$mean[,"estimates"],
                row.names(fit.2ht.5k$parameters$mean))

# for (file in load_files) {
#   load(file = file)
# }
15*(1-par["do"])*par["g5"]*par["s_g1"]*par["s_g2"]*par["s_g3"]*par["s_g4"]



sel2 <- select.mpt(list(fit.2ht.5k, fit.sdt.5k, fit.aj.5k.3clk, fit.aj.5k.3gk))

sel <- transform(sel2,
                 "p%smaller.05" = (p.smaller.05 / 47) * 100,
                 "AIC%best" = (AIC.best / 47) * 100)[,
                                                     c("model", "n.parameters",
                                                       "p.smaller.05", "p%smaller.05",
                                                       "AIC.best", "AIC%best")]
write_xlsx(sel, "sel_output_5k.xlsx")

# Funciones de formato
format_p <- function(p_values) {
  sapply(p_values, function(p) {
    if (p < 0.001) {
      "<.001"
    } else {
      sub("^0", "", format(round(p, 3), nsmall = 3)) # Sin 0 inicial y 3 decimales
    }
  })
}

format_x2_aic <- function(values) {
  format(round(values, 2), nsmall = 2) # Redondeo a 2 decimales
}

# Creación de la tabla con nuevos modelos
tabla2 <- data.frame(   # Table B.2
  Id = 1:47,
  `X2(df=51)_2HT` = format_x2_aic(fit.2ht.5k$goodness.of.fit$individual$G.Squared),
  p_2HT = format_p(fit.2ht.5k$goodness.of.fit$individual$p.value),
  AIC_2HT = format_x2_aic(fit.2ht.5k$information.criteria$individual$AIC),
  
  `X2(df=57)_SDT` = format_x2_aic(fit.sdt.5k$goodness.of.fit$individual$G.Squared),
  p_SDT = format_p(fit.sdt.5k$goodness.of.fit$individual$p.value),
  AIC_SDT = format_x2_aic(fit.sdt.5k$information.criteria$individual$AIC),
  
  `X2(df=50)_Aj5k.3clk` = format_x2_aic(fit.aj.5k.3clk$goodness.of.fit$individual$G.Squared),
  p_AJcj = format_p(fit.aj.5k.3clk$goodness.of.fit$individual$p.value),
  AIC_AJcj = format_x2_aic(fit.aj.5k.3clk$information.criteria$individual$AIC),
  
  `X2(df=50)_Aj5k.3gk` = format_x2_aic(fit.aj.5k.3gk$goodness.of.fit$individual$G.Squared),
  p_AJgk = format_p(fit.aj.5k.3gk$goodness.of.fit$individual$p.value),
  AIC_AJgk = format_x2_aic(fit.aj.5k.3gk$information.criteria$individual$AIC)
)

# Guardar la tabla en un archivo Excel
write_xlsx(tabla2, "tabla1_5k.xlsx")


#-------------------------------------------------------------------------------
# Goodness of fit
#-------------------------------------------------------------------------------

# Goodness of fit for 5j 2CL models
Gof.5j.2CL <- data.frame(
  "AJcj" = cbind.data.frame(
    fit.aj.5cj.2CL$goodness.of.fit$individual,
    aic = fit.aj.5cj.2CL$information.criteria$individual$AIC
  ),
  "AJgj" = cbind.data.frame(
    fit.aj.5gj.2CL$goodness.of.fit$individual,
    aic = fit.aj.5gj.2CL$information.criteria$individual$AIC
  ),
  "2HT" = cbind.data.frame(
    fit.2ht.5j.3CL$goodness.of.fit$individual,
    aic = fit.aj.5cj.2CL$information.criteria$individual$AIC
  ),
  "SDT" = cbind.data.frame(
    fit.sdt.5j.2CL$goodness.of.fit$individual,
    aic = fit.aj.5cj.2CL$information.criteria$individual$AIC
  )
)

# Goodness of fit for 3j 2CL models
Gof.3j.2CL <- data.frame(
  "AJcj" = cbind.data.frame(
    fit.aj.3clk$goodness.of.fit$individual,
    aic = fit.aj.3clk$information.criteria$individual$AIC
  ),
  "AJgj" = cbind.data.frame(
    fit.aj.3gk$goodness.of.fit$individual,
    aic = fit.aj.3gk$information.criteria$individual$AIC
  ),
  "2HT" = cbind.data.frame(
    fit.2ht$goodness.of.fit$individual,
    aic = fit.2ht$information.criteria$individual$AIC
  ),
  "SDT" = cbind.data.frame(
    fit.sdt$goodness.of.fit$individual,
    aic = fit.sdt$information.criteria$individual$AIC
  )
)

# Save goodness of fit dataframes to Excel
write_xlsx(x = Gof.5j.2CL, path = "Gof.5j.2CL.xlsx")
write_xlsx(x = Gof.3j.2CL, path = "Gof.3j.2CL.xlsx")


#-------------------------------------------------------------------------------
# 6.4. Cross Validation
#-------------------------------------------------------------------------------

datos <- read.delim("data_Juola.txt")
datos$Confidence_level
datos <- datos[(datos$Target_frequency != 1 & datos$Target_frequency != 5),]
datos[(datos$Confidence_level==3),"Confidence_level"] <- 2 
datos$Target_frequency <- datos$Target_frequency - 1



# Initial time
ini_time <- Sys.time()

#-----
# 3.a. Prepare for parallel
#-----

n.cores <- parallel::detectCores(logical = FALSE)-1
cl <- parallel::makeCluster(n.cores, type = "SOCK") 
doParallel::registerDoParallel(cl, cores = n.cores)


pm <- foreach::foreach(j = 1:k, # performance metrics (pm)
                       .packages = c("caret", "dplyr", "tidyr", "MPTinR"), 
                       .combine = rbind, 
                       .inorder = FALSE
) %dopar% {
  try({ 
    
    #-----
    # 3.1. Establish randomization seed
    #-----
    
    set.seed(1320 + j) # Unique and identifiable seed for each condition
    
    #-----
    # 3.2. Specify condition
    
    datos_id <- datos[datos$Id == 1, ]
    indices <- createDataPartition(y = datos_id$Category, p = 0.5, list = FALSE)
    train_set <- datos_id[indices, ]
    test_set  <- datos_id[-indices, ]
    
    for (i in 2:max(datos$Id)) {
      
      datos_id <- datos[datos$Id == i, ]
      indices <- createDataPartition(y = datos_id$Category, p = 0.5, list = FALSE)
      train_set2 <- datos_id[indices, ]
      test_set2  <- datos_id[-indices, ]
      train_set <- rbind.data.frame(train_set, train_set2)
      test_set <- rbind.data.frame(test_set, test_set2)
      
    }
    
    train_set <- change_to_freq(train_set)
    test_set <- change_to_freq(test_set)
    
    train.2ht<-  fit.model(model.filename = textConnection(m.2ht), data = train_set, 
                           multicore = "individual", nCPU = 2, sfInit = TRUE,
                           n.optim =  20,
                           ci = 95, 
                           use.gradient = F)
    
    train.sdt<-  fit.model(model.filename = textConnection(m.sdt), data = train_set, 
                         multicore = "individual", nCPU = 2, sfInit = TRUE,
                         n.optim =  20,
                         ci = 95, 
                         lower.bound= c(rep(-Inf,3), 0,0, 0.1, rep(0,4), 0.1),
                         upper.bound= c(rep(Inf, 6), rep(1, 4), Inf),
                         use.gradient = F)
    
    train.aj <-   fit.model(model.filename = textConnection(m.aj), data = train_set, 
                          multicore = "individual", nCPU = 2, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          output = "full",
                          lower.bound= c(0.1, rep(-Inf,3), 0.1, rep(0,14), 0.1),
                          upper.bound= c(rep(Inf,5), rep(1,14), Inf))
    
    test.2ht<-  fit.model(model.filename = textConnection(m.2ht), data = test_set, 
                          multicore = "individual", nCPU = 2, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          use.gradient = F)
    test.aj <-   fit.model(model.filename = textConnection(m.aj), data = test_set, 
                           multicore = "individual", nCPU = 2, sfInit = TRUE,
                           n.optim =  20,
                           ci = 95, 
                           output = "full",
                           lower.bound= c(0.1, rep(-Inf,3), 0.1, rep(0,14), 0.1),
                           upper.bound= c(rep(Inf,5), rep(1,14), Inf))
    
    test.sdt<-  fit.model(model.filename = textConnection(m.sdt), data = test_set, 
                          multicore = "individual", nCPU = 2, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          lower.bound= c(rep(-Inf,3), 0,0, 0.1, rep(0,4), 0.1),
                          upper.bound= c(rep(Inf, 6), rep(1, 4), Inf),
                          use.gradient = F)
    
    
    
    res <- c(select.mpt(list(test.2ht, test.sdt, test.aj))[, "AIC.best"], 
                  select.mpt(list(train.2ht, train.sdt, train.aj))[, "AIC.best"])
    names(res) <- c("test.2ht", "test.sdt", "test.aj",
                        "train.2ht", "train.sdt", "train.aj")
    res # Final output tibbles
  })
}

fin_time <- Sys.time(); 
print(elapse_time <- fin_time - ini_time)

parallel::stopCluster(cl)



library(dplyr)
library(tidyr)
library(writexl)

resumen_datos <- datos %>%
  group_by(id, x) %>%
  summarize(f_cat = n(),           
            m_y = mean(y, na.rm = TRUE))  


tabla_final <- resumen_datos %>%
  pivot_wider(
    names_from = x,               
    values_from = c(f_cat, m_y),  
    names_glue = "{.value}_x{x}"  
  )


table(datos$Id, datos$Confidence_level, datos$Confidence_level)
writexl::write_xlsx(tabla_final, "tabla_final.xlsx")
tabla_cnf <- (table(datos$id, datos$x_cnf))

df_tabla_cnf <- as.data.frame.matrix(tabla_cnf)


writexl::write_xlsx(df_tabla_cnf, "df_tabla_cnf.xlsx")
writexl::write_xlsx(as.data.frame.matrix(par.3clk), "par.3clk.xlsx")



#-------------------------------------------------------------------------------
# 7. Data Plots
#-------------------------------------------------------------------------------

#-----
# 7.1. Observed Proportion
#-----
library(dplyr)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)

long_counts <- as.data.frame(d.data.cl.rt) %>%
  mutate(subject = row_number()) %>%                           
  pivot_longer(
    cols = -subject,
    names_to = "var",
    values_to = "count"
  ) %>%
  # 
  tidyr::extract(
    var,
    into  = c("response","k","conf","speed"),
    regex = "^([a-z]+)-k([123])-(high|low)([12])$"
  ) %>%
  mutate(
    target_freq = recode(k, `1` = "65%", `2` = "50%", `3` = "35%"),
    stimulus_type = if_else(response %in% c("hit","miss"), "Target", "Lure"),
    confidence = if_else(conf == "high", "High", "Low"),
    rt_bin     = if_else(speed == "1", "Fast", "Slow")
  )


totals_by_subject <- long_counts %>%
  group_by(subject, k, target_freq, stimulus_type) %>%
  summarise(total_stim = sum(count), .groups = "drop")


props_by_subject <- long_counts %>%
  left_join(totals_by_subject,
            by = c("subject","k","target_freq","stimulus_type")) %>%
  mutate(prop = if_else(total_stim > 0, count / total_stim, NA_real_))


summary_prop <- props_by_subject %>%
  group_by(response, target_freq, stimulus_type, confidence, rt_bin) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop   = sd(prop,  na.rm = TRUE),
    n         = sum(!is.na(prop)),
    se_prop   = sd_prop / sqrt(n),
    ci95_prop = 1.96 * se_prop,
    .groups = "drop"
  ) %>%
  mutate(
    response = factor(response,
                      levels = c("hit", "miss", "fa", "cr"),
                      labels = c("Hit", "Miss", "False Alarm", "Correct Rejection")),
    target_freq = factor(target_freq,
                         levels = c("65%", "50%", "35%"),
                         labels = c("65% Target", "50% Target", "35% Target"))
  )

p <- ggplot(summary_prop,
            aes(x = interaction(confidence, rt_bin),
                y = mean_prop, fill = confidence)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean_prop - ci95_prop,
                    ymax = mean_prop + ci95_prop),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  facet_grid(target_freq ~ response) +
  scale_fill_grey(start = 0.3, end = 0.7, name = "Confidence") +
  labs(
    title = NULL,
    x = "CL × RT",
    y = "Mean Proportion (95% CI)"
  ) +
  theme_bw(base_size = 10, base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

# ggsave("mean_proportion_barplot.png", plot = p, width = 8, height = 5, dpi = 300)



summary_hits_cr_highconf <- summary_prop %>%
  filter(response %in% c("Hit", "Correct Rejection"),
         confidence == "High") %>%
  mutate(
    ci_lower = mean_prop - 1.96 * se_prop,
    ci_upper = mean_prop + 1.96 * se_prop
  )

# CI
p2 <- ggplot(summary_hits_cr_highconf,
             aes(x = target_freq,
                 y = mean_prop,
                 group = rt_bin,
                 shape = rt_bin,
                 linetype = rt_bin)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.1, linewidth = 0.3) +
  geom_line() +
  facet_wrap(~ response, nrow = 1) +
  scale_y_continuous(name = "Mean Proportion", limits = c(0.05, 0.45)) +
  scale_x_discrete(name = "Target Frequency Condition") +
  scale_shape_manual(values = c("Fast" = 16, "Slow" = 1)) +
  scale_linetype_manual(values = c("Fast" = "solid", "Slow" = "dashed")) +
  theme_bw(base_size = 8, base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = c(0.98, 0.65),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10),
    strip.text.x = element_text(size = 10)
  ) +
  labs(
    shape = "Response Time",
    linetype = "Response Time"
  )

p2


# ggsave("selected_mean_proportion_barplot.png", plot = p2, width = 8, height = 3, dpi = 300)



#-----
# 7.1. Predicted vs Observed Proportion
#-----

# path: R/plot_obs_pred_high_function.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

plot_obs_pred_high <- function(fit, y_limits = c(0.05, 0.45)) {
  # Why: ensure required slots exist and are matrices
  if (is.null(fit$data$observed$individual) || is.null(fit$data$predicted$individual)) {
    stop("fit must contain $data$observed$individual and $data$predicted$individual.")
  }
  m_obs  <- fit$data$observed$individual
  m_pred <- fit$data$predicted$individual
  if (!is.matrix(m_obs) || !is.matrix(m_pred)) stop("observed/predicted must be matrices.")
  
  # Why: unify predicted colnames to observed pattern when shape matches
  if (!identical(colnames(m_pred), colnames(m_obs)) && ncol(m_pred) == ncol(m_obs)) {
    colnames(m_pred) <- colnames(m_obs)
  }
  
  # Keep only Hit/CR high1/high2 for k1..k3
  keep_pat <- "^(hit|cr)-k[123]-high[12]$"
  obs_high  <- m_obs [, grepl(keep_pat, colnames(m_obs )), drop = FALSE]
  pred_high <- m_pred[, grepl(keep_pat, colnames(m_pred)), drop = FALSE]
  if (ncol(obs_high) == 0 || ncol(pred_high) == 0) {
    stop("No columns matched '^(hit|cr)-k[123]-high[12]$'.")
  }
  
  # Denominators per column
  make_denoms <- function(nms) {
    k_num  <- as.integer(sub(".*-k([123])-.*", "\\1", nms))
    is_hit <- grepl("^hit", nms)
    den_hit <- c(65, 50, 35) # k1,k2,k3 targets
    den_cr  <- c(35, 50, 65) # k1,k2,k3 lures
    ifelse(is_hit, den_hit[k_num], den_cr[k_num])
  }
  den_obs  <- make_denoms(colnames(obs_high))
  den_pred <- make_denoms(colnames(pred_high))
  
  # Proportions per subject/column
  obs_prop  <- sweep(obs_high,  2, den_obs,  "/")
  pred_prop <- sweep(pred_high, 2, den_pred, "/")
  
  # Long format with labels
  to_long <- function(prop_mat, source_label) {
    as.data.frame(prop_mat) |>
      mutate(Subject = dplyr::row_number()) |>
      tidyr::pivot_longer(-Subject, names_to = "col", values_to = "prop") |>
      mutate(
        Response = if_else(str_starts(col, "hit"), "Hit", "Correct Rejection"),
        K = factor(stringr::str_match(col, "-k([123])-")[, 2],
                   levels = c("1","2","3"),
                   labels = c("65% Target","50% Target","35% Target")),
        Speed = if_else(str_detect(col, "high1"), "Fast", "Slow"),
        Source = source_label
      ) |>
      select(Subject, Response, K, Speed, Source, prop)
  }
  long_obs  <- to_long(obs_prop,  "Observed")
  long_pred <- to_long(pred_prop, "Predicted")
  long_all  <- bind_rows(long_obs, long_pred)
  
  # Summary stats
  summary_both <- long_all |>
    group_by(Response, K, Speed, Source) |>
    summarise(
      n = sum(!is.na(prop)),
      mean_prop = mean(prop, na.rm = TRUE),
      sd_prop   = sd(prop,   na.rm = TRUE),
      se_prop   = sd_prop / sqrt(n),
      ci_lower  = mean_prop - 1.96 * se_prop,
      ci_upper  = mean_prop + 1.96 * se_prop,
      .groups = "drop"
    ) |>
    mutate(
      # Why: enforce panel and legend order
      Response = factor(Response, levels = c("Hit","Correct Rejection")),
      Speed    = factor(Speed, levels = c("Fast","Slow")),
      Source   = factor(Source, levels = c("Observed","Predicted"))
    )
  
  # Plot
  p_both <- ggplot(
    summary_both,
    aes(x = K, y = mean_prop,
        group = interaction(Source, Speed),
        color = Source, linetype = Source, shape = Speed)
  ) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.10, linewidth = 0.3) +
    geom_line() +
    geom_point(size = 3) +
    facet_wrap(~ Response, nrow = 1) +
    scale_y_continuous(name = "Mean Proportion", limits = y_limits) +
    scale_x_discrete(name = "Target Frequency Condition") +
    scale_shape_manual(values = c("Fast" = 16, "Slow" = 1)) +
    theme_bw(base_size = 8, base_family = "serif") +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 10)
    ) +
    labs(
      color = "Source",
      linetype = "Source",
      shape = "Response Time"
    )
  
  list(summary = summary_both, plot = p_both)
}

# --- Example usage ---
res_aj <- plot_obs_pred_high(fit.aj.3clk)
print(res_aj$plot)
ggsave("observed_vs_predicted_AJ.png", res_aj$plot, width = 8, height = 3, dpi = 300)

res_2ht <- plot_obs_pred_high(fit.2ht)
print(res_2ht$plot)
ggsave("observed_vs_predicted_2HT.png", res_2ht$plot, width = 8, height = 3, dpi = 300)

res_sdt <- plot_obs_pred_high(fit.sdt)
print(res_sdt$plot)
ggsave("observed_vs_predicted_SDT.png", res_sdt$plot, width = 8, height = 3, dpi = 300)


#-----
# 7.3. Full
#-----

wd <- getwd()
df <- read.table( paste0(wd, "/data/data_juola.txt"), header = TRUE, sep = "\t")


df_clean <- df %>%
  filter(Target_frequency %in% c(2, 3, 4)) %>%
  mutate(
    trial_type = case_when(
      Type_of_trial == 1 ~ "Hit",
      Type_of_trial == 2 ~ "Miss",
      Type_of_trial == 3 ~ "FA",
      Type_of_trial == 4 ~ "CR"
    ),
    confidence_bin = case_when(
      Confidence_level == 1 ~ "High",
      Confidence_level %in% c(2, 3) ~ "Low"
    ),
    target_freq_label = case_when(
      Target_frequency == 2 ~ "65%",
      Target_frequency == 3 ~ "50%",
      Target_frequency == 4 ~ "35%"
    )
  )

mean_freq_per_subject <- df_clean %>%
  group_by(Id, trial_type, confidence_bin, target_freq_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(trial_type, confidence_bin, target_freq_label) %>%
  summarise(mean_n = mean(n), .groups = "drop")

# Paso 3: calcular RTs como antes
summary_rt <- df_clean %>%
  group_by(trial_type, confidence_bin, target_freq_label) %>%
  summarise(
    median_rt = median(RT, na.rm = TRUE),
    ric_rt = IQR(RT, na.rm = TRUE),
    .groups = "drop"
  )

summary_final <- left_join(mean_freq_per_subject, summary_rt,
                           by = c("trial_type", "confidence_bin", "target_freq_label"))

summary_table_ordered <- summary_final %>%
  arrange(target_freq_label, trial_type, confidence_bin) %>%
  mutate(
    mean_n = round(mean_n, 1),
    median_rt = round(median_rt),
    ric_rt = round(ric_rt),
    target_freq_label = factor(target_freq_label, levels = c("35%", "50%", "75%"))
  ) %>%
  select(target_freq_label, trial_type, confidence_bin, mean_n, median_rt, ric_rt)

library(dplyr)
library(knitr)
library(kableExtra)

# summary_table_ordered %>%
#   select(-target_freq_label) %>%
#   kable(
#     format = "html",
#     col.names = c("Response Category", "Confidence", "Mean Freq.", "Median RT (ms)", "IQR (ms)"),
#     align = "lcccc"
#   ) %>%
#   kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE, position = "center") %>%
#   row_spec(0, bold = TRUE) %>%
#   group_rows("Target Freq.: 35%", 1, 8) %>%
#   group_rows("Target Freq.: 50%", 9, 16) %>%
#   group_rows("Target Freq.: 75%", 17, 24)
# 
# library(webshot2)
# library(magick)
# install.packages("magick")
# library(kableExtra)
# 
# save_kable(summary_table_ordered, file = "descript_freq_medianRT.png")

descript_freq_medianRT <- summary_table_ordered %>%
  select(-target_freq_label) %>%
  kable(
    format = "html",
    col.names = c("Response Category", "Confidence", "Mean Freq.", "Median RT (ms)", "IQR (ms)"),
    align = "lcccc"
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center") %>%
  row_spec(0, bold = TRUE) %>%
  group_rows("Target Freq.: 35%", 1, 8) %>%
  group_rows("Target Freq.: 50%", 9, 16) %>%
  group_rows("Target Freq.: 75%", 17, 24)

# Guardar como imagen
save_kable(descript_freq_medianRT, file = "descript_freq_medianRT.png")

-----
  # 7.2. Table: Median RT and IQR by Trial Type, Confidence, and Target Frequency
  #-----
wd <- getwd()
df <- read.table(paste0(wd, "/data/data_juola.txt"), header = TRUE, sep = "\t")

# 2. Cargar paquetes
library(dplyr)
library(ggplot2)

# # 3. Preparar datos para el gráfico
# plot_data <- df %>%
#   filter(Target_frequency %in% c(2, 3, 4)) %>%
#   mutate(
#     trial_type = case_when(
#       Type_of_trial == 1 ~ "Hit",
#       Type_of_trial == 2 ~ "Miss",
#       Type_of_trial == 3 ~ "FA",
#       Type_of_trial == 4 ~ "CR"
#     ),
#     confidence_bin = case_when(
#       Confidence_level == 1 ~ "High",
#       Confidence_level %in% c(2, 3) ~ "Low"
#     ),
#     target_freq_label = case_when(
#       Target_frequency == 2 ~ "65%",
#       Target_frequency == 3 ~ "50%",
#       Target_frequency == 4 ~ "35%"
#     ),
#     condition = paste0(tolower(trial_type), "_", tolower(confidence_bin)),
#     group = ifelse(trial_type %in% c("CR", "FA"), "Noise", "Signal")
#   ) %>%
#   group_by(trial_type, confidence_bin, target_freq_label, condition, group) %>%
#   summarise(
#     median_rt = median(RT, na.rm = TRUE),
#     ric_rt = IQR(RT, na.rm = TRUE),
#     .groups = "drop"
#   )

# 4. Etiquetas legibles para condiciones
condition_order <- c(
  "cr_high", "cr_low", "fa_low", "fa_high",
  "miss_high", "miss_low", "hit_low", "hit_high"
)
condition_labels <- c(
  "cr_high" = "CR High",
  "cr_low" = "CR Low",
  "fa_low" = "FA Low",
  "fa_high" = "FA High",
  "miss_high" = "Miss High",
  "miss_low" = "Miss Low",
  "hit_low" = "Hit Low",
  "hit_high" = "Hit High"
)

plot_data <- plot_data %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    condition_pretty = factor(condition, labels = condition_labels),
    target_freq_label = factor(target_freq_label, levels = c("35%", "50%", "75%"))
  )

# 5. Etiquetas para facet rows
facet_labels <- c(
  "35%" = "Target Frequency: 35%",
  "50%" = "Target Frequency: 50%",
  "75%" = "Target Frequency: 75%"
)

# 6. Generar el gráfico
ggplot(plot_data, aes(x = condition_pretty, y = median_rt, color = group, group = group)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = median_rt - ric_rt/2, ymax = median_rt + ric_rt/2), width = 0.2) +
  facet_grid(rows = vars(target_freq_label), labeller = as_labeller(facet_labels)) +
  scale_color_manual(values = c("Noise" = "black", "Signal" = "gray50")) +
  labs(
    x = "Category x CL",
    y = "Median RT (ms)",
    color = "Stimuly"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

table(plot_data$condition_pretty, plot_data$group)
df %>% count(Type_of_trial, Confidence_level, Stimulus, Response)
