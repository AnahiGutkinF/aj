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
              "data.table", "parallel", "writexl", "caret", "snowfall", "dplyr", 
              "knitr", "kableExtra", "stringr")

lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

source("aj.fun.R")

data_dir <- path.expand("~/GitHub/aj/data")
nml_dir <- path.expand("~/GitHub/aj/nml")
fit_dir <- path.expand("~/GitHub/aj/main_fits")
fig_dir <- path.expand("~/GitHub/aj/figures")

# -------------------------------------------------------------------------------
# 2. Conditions & Data
#-------------------------------------------------------------------------------

objetos <- ls()
cond <- c("aj.3clk","aj.3gk", "sdt", "2ht")
f <- paste0("fit.", cond) 
m <- paste0("m.", cond) 

conditions <- cbind("fit" = f, "model"=m); rm(cond, f, m, objetos)
load(file.path( data_dir, "d.data.cl.rt.RData")) #load RT binned data 
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
                     use.gradient = F)
save(fit.2ht, file = file.path(fit_dir, "fit.2ht.RData"))

fit.sdt<-  fit.model(model.filename = m.sdt, data = d.data.cl.rt, 
                     multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                     n.optim =  20,
                     ci = 95, 
                     lower.bound= c(-Inf, 0, 0, 0.1, rep(0,6), 0.1),
                     upper.bound= c(rep(Inf, 4), rep(1, 2), rep(Inf,2), rep(1,2), Inf),
                     use.gradient = F)
save(fit.sdt, file = file.path(fit_dir, "fit.sdt.RData"))


fit.aj <-   fit.model(model.filename = m.aj, data = d.data.cl.rt, 
                      multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                      n.optim =  20,
                      ci = 95, 
                      output = "full",
                      lower.bound= c(0.1, rep(-Inf,1), 0.1, rep(0,16), 0.1),
                      upper.bound= c(rep(Inf,3), rep(1,5), Inf, Inf, rep(1,9), Inf)
                      )
save(fit.aj, file = file.path(fit_dir, "fit.aj.RData"))

fit.aj.3clk <-  fit.model(model.filename = m.aj, data = d.data.cl.rt, 
                          multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          output = "full",
                          lower.bound= c(0.1, rep(-8, 1), 0.1, rep(0,14), 0.1),
                          upper.bound= c(rep(Inf,3), rep(1, 5), rep(Inf,2), rep(1,7), 1),
                          restrictions.filename = list("s_g1=1", "s_g2=1")
                          )
save(fit.aj.3clk, file = file.path(fit_dir, "fit.aj.3clk.RData"))

# check.mpt(m.aj, restrictions.filename = list("s_cl2=1", "s_cl3=1"))
fit.aj.3gk <-  fit.model(model.filename = m.aj.3gk, data = d.data.cl.rt, 
                         multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                          n.optim =  20,
                          ci = 95, 
                          output = "full",
                          lower.bound= c(0.1, -8, 0.1, rep(0,14), 0.1),
                          upper.bound= c(rep(Inf, 3), rep(1, 14), Inf),
                          restrictions.filename = list("s_cl2=1", "s_cl3=1"),
                          use.gradient = T)
save(fit.aj.3gk, file.path(fit_dir, file = "fit.aj.3gk.RData"))


# Loas main fits
load_files <- c("fit.2ht.RData", "fit.sdt.RData", "fit.aj.RData", 
                "fit.aj.3clk.RData", "fit.aj.3gk.RData")

for (file in load_files) {
  load(file = file.path(fit_dir, file))
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

tabla1 <- data.frame(   # Table B.3
  Id = 1:n_sub,
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

writexl::write_xlsx(tabla1, path = file.path(fig_dir, "tableB3_3k.xlsx"))

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

save(fit.aj.3clk.r1, file = file.path(fit_dir, "fit.aj.3clk.r1.RData"))
load(file.path(fit_dir, "fit.aj.3clk.r1.RData"))

nested.cj <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.r1)
nested.gj <- test.nested(model.gen = fit.aj.3gk , model.rest = fit.aj.3clk.r1)

sum(nested.cj$p.val<0.05)/n_sub
sum(nested.gj$p.val<0.05)/n_sub

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

table_nested1 <- data.frame( # Table D.1
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

save(fit.aj.3clk.rL, file = file.path(fit_dir,"fit.aj.3clk.rL.RData"))
load(file.path(fit_dir,"fit.aj.3clk.rL.RData"))

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

save(fit.aj.3g.rL, file = file.path(fit_dir,"fit.aj.3g.rL.RData"))
load(file.path(fit_dir,"fit.aj.3g.rL.RData"))

nested.cj.L <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.rL)
nested.gj.L <- test.nested(model.gen = fit.aj.3gk , model.rest = fit.aj.3g.rL)

sum(nested.cj.L$p.val<.05)/n_sub
sum(nested.gj.L$p.val<.05)/n_sub

table_nested2 <- data.frame( # Table D.2
  Id = 1:47,
  
  `G2 (df=3)_AJcj` = nested.cj.L$G2,
  `p_AJcj` = format_p(nested.cj.L$p.val),
  
  `G2 (df=3)_AJgj` = nested.gj.L$G2,
  `p_AJgj` = format_p(nested.gj.L$p.val)
)

#-------------------------------------------------------------------------------
# 3.2.4 Testing A-J assumptions about CLs (H parameters)
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
save(fit.aj.3clk.rH, file = file.path(fit_dir,"fit.aj.3clk.rH.RData"))
load(file.path(fit_dir,"fit.aj.3clk.rH.RData"))

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
save(fit.aj.3gk.rH, file = file.path(fit_dir,"fit.aj.3gk.rH.RData"));
load(file.path(fit_dir,"fit.aj.3gk.rH.RData"))

nested.cj.H <- test.nested(model.gen = fit.aj.3clk, model.rest = fit.aj.3clk.rH)
nested.gj.H <- test.nested(model.gen = fit.aj.3gk, model.rest = fit.aj.3gk.rH )

sum(nested.cj.H$p.val<0.05)/n_sub
sum(nested.gj.H$p.val<0.05)/n_sub
    
    table_nested3 <- data.frame( # Table D.3
      Id = 1:n_sub,
      
      `G2 (df=4)_AJcj` = nested.cj.H$G2,
      `p_AJcj` = format_p(nested.cj.H$p.val),
      
      `G2 (df=4)_AJgj` = nested.gj.H$G2,
      `p_AJgj` = format_p(nested.gj.H$p.val)
    )

#-------------------------------------------------------------------------------
# 3.3. p-values plots (Figure 7 - 9)
#-------------------------------------------------------------------------------

  plot_aj_logp_histogram <- function(gj_obj, cj_obj, label = "AJ") {
    name_g <- paste0(label, "g")
    name_c <- paste0(label, "c")
    
    gj_vals <- as.numeric(unlist(gj_obj$p.val))
    cj_vals <- as.numeric(unlist(cj_obj$p.val))
    
    gj_vals[gj_vals == 0] <- 0.001
    cj_vals[cj_vals == 0] <- 0.001
    
    gj_df <- data.frame(
      log_pval = log(gj_vals),
      test = name_g
    )
    cj_df <- data.frame(
      log_pval = log(cj_vals),
      test = name_c
    )
    
    plot_df <- dplyr::bind_rows(gj_df, cj_df)
    threshold_log <- log(0.05)
    
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
        legend.position = "none"  # hide legend
      )
  }
  
  p_AJ_hist    <- plot_aj_logp_histogram(nested.gj, nested.cj, label = "AJ") # Fig7
  p_AJ_hist_L  <- plot_aj_logp_histogram(nested.gj.L, nested.cj.L, label = "AJ") # Fig8
  p_AJ_hist_H  <- plot_aj_logp_histogram(nested.gj.H, nested.cj.H, label = "AJ") # Fig9
  
  ggsave(file.path(fig_dir,"aj_hist_hypothesis_testing.png"),
         plot = p_AJ_hist,
         width = 8, height = 4, dpi = 300, units = "in")
  
  ggsave(file.path(fig_dir,"aj__hist_hypothesis_testing_H.png"),
         plot = p_AJ_hist_H,
         width = 8, height = 4, dpi = 300, units = "in")
  
  ggsave(file.path(fig_dir,"aj__hist_hypothesis_testing_L.png"),
         plot = p_AJ_hist_L,
         width = 8, height = 4, dpi = 300, units = "in")
  
#-------------------------------------------------------------------------------
# 4. Import fits 
#-------------------------------------------------------------------------------

archivos_fit <- list.files(fit_dir, pattern = "^fit\\.")
for (archivo in archivos_fit) {
  ruta_archivo <- file.path(fit_dir, archivo)
  load(ruta_archivo, .GlobalEnv)  # Carga el objeto en el entorno global
  cat("Imported data:", archivo, "\n")
  }

#import nml penality

load(file.path(nml_dir, "nmlaj.RData"))
load(file.path(nml_dir, "nml2ht.RData"))
load(file.path(nml_dir, "nmlsdt.RData"))

#-------------------------------------------------------------------------------
# 5. Grouping par estimations & Mean par
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
# 6. Restructure data with 5 target conditions 
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
d.juola <- read.delim(file.path(data_dir, "data_Juola.txt"))
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
save(data.5j.3CL, file = file.path(data_dir, "data.5j.3CL.RData"))
save(data.5j.2CL, file = file.path(data_dir, "data.5j.2CL.RData"))

#-------------------------------------------------------------------------------
# 7. Model fitting for 5 %target conditions
#-------------------------------------------------------------------------------

# check.mpt(m.aj.5k, restrictions.filename = list("s_g1=1", "s_g2=1",
#                                                 "s_g3=1", "s_g4=1"))
# check.mpt(m.aj.5k, restrictions.filename = list("s_cl2=0", "s_cl3=0",
#                                                 "s_cl4=0", "s_cl5=0"))

n.cores <- parallel::detectCores()-1

fit.2ht.5k <- fit.model(model.filename = m.2ht.5k, data = data.5j.2CL, 
                         multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                         n.optim =  20,
                         ci = 95, 
                         use.gradient = F,
                         )
save(fit.2ht.5k, file = file.path(fit_dir, "fit.2ht.5k.RData"))

fit.aj.5k.3clk <- fit.model(model.filename = m.aj.5k, data = data.5j.2CL, 
                           multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                           n.optim =  20,
                           ci = 95,  
                           lower.bound = c(0.1, -8, 0.1, rep(0, 16), 0.1),
                           upper.bound = c(rep(Inf, 3), rep(1, 5), rep(1,4), rep(1,7), Inf),
                           use.gradient = F,
                           restrictions.filename = list("s_g1=1", "s_g2=1",
                                                        "s_g3=1", "s_g4=1"),
                           )
save(fit.aj.5k.3clk, file = file.path(fit_dir, "fit.aj.5k.3clk.RData"))

fit.aj.5k.3gk <- fit.model(model.filename = m.aj.5k, data = data.5j.2CL, 
                       multicore = "individual", nCPU = n.cores, sfInit = TRUE,
                       n.optim =  20,
                       ci = 95,  
                       lower.bound = c(0.1, -8, 0.1, rep(0, 16), 0.1),
                       upper.bound = c(rep(Inf, 3), rep(1, 16), Inf),
                       use.gradient = F,
                       restrictions.filename = list("s_cl2=0", "s_cl3=0",
                                                    "s_cl4=0", "s_cl5=0"),
                       )
save(fit.aj.5k.3gk, file = file.path(fit_dir, "fit.aj.5k.3gk.RData")) 

fit.sdt.5k <- fit.model(model.filename = m.sdt.5k, data = data.5j.2CL, 
                            multicore = "n.optim", sfInit = TRUE, ci = 95, 
                            lower.bound = c(-8, rep(0,2), 0.1, rep(0, 8), 0.1),
                            upper.bound = c(rep(Inf, 4), rep(1, 2), rep(Inf,4), rep(1,2), Inf),
                            use.gradient = FALSE
                            )
save(fit.sdt.5k, file = file.path(fit_dir, "fit.sdt.5k.RData"))


# Load 5 categ. fits
load_files <- c("fit.2ht.5k.RData", "fit.sdt.5k.RData",
                "fit.aj.5k.3gk.RData", "fit.aj.5k.3clk.RData")

for (file in load_files) {
  load(file = file.path(fit_dir, file))
}

#Appendix B (Table B.1)
par <- setNames(fit.2ht.5k$parameters$mean[,"estimates"],
                row.names(fit.2ht.5k$parameters$mean))


# Expected Freq for “guessing–high-confidence–fast
# 15*(1-par["do"])*par["g5"]*par["s_g1"]*par["s_g2"]*par["s_g3"]*par["s_g4"]

sel2 <- select.mpt(list(fit.2ht.5k, fit.sdt.5k, fit.aj.5k.3clk, fit.aj.5k.3gk))

sel <- transform(sel2,
                 "p%smaller.05" = (p.smaller.05 / n_sub) * 100,
                 "AIC%best" = (AIC.best / n_sub) * 100)[,
                                                     c("model", "n.parameters",
                                                       "p.smaller.05", "p%smaller.05",
                                                       "AIC.best", "AIC%best")]
# Format Table B.2
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

# Built Table B.2
tabla2 <- data.frame(   # Table B.2
  Id = 1:n_sub,
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
write_xlsx(tabla2, path = file.path(fig_dir , "tablaB2_5k.xlsx"))


#-------------------------------------------------------------------------------
# 8. Data Plots
#-------------------------------------------------------------------------------

#-----
# 8.1. Observed Proportion (Figure 6)
#-----

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

print(p) #Fig6
ggsave( file.path(fig_dir, "mean_proportion_barplot.png"),
        plot = p, width = 8, height = 5, dpi = 300)

#-----
# 8.2. Predicted vs Observed Proportion (Fig 6A-6C)
#-----

plot_obs_pred_high <- function(fit, y_limits = c(0.05, 0.5)) {
   # Ensure requirements
  if (is.null(fit$data$observed$individual) || is.null(fit$data$predicted$individual)) {
    stop("fit must contain $data$observed$individual and $data$predicted$individual.")
  }
  m_obs  <- fit$data$observed$individual
  m_pred <- fit$data$predicted$individual
  if (!is.matrix(m_obs) || !is.matrix(m_pred)) stop("observed/predicted must be matrices.")
  
  if (!identical(colnames(m_pred), colnames(m_obs)) && ncol(m_pred) == ncol(m_obs)) {
    colnames(m_pred) <- colnames(m_obs)
  }
  
  # Keep only Correct High confidence responses 
  keep_pat <- "^(hit|cr)-k[123]-high[12]$"
  obs_high  <- m_obs [, grepl(keep_pat, colnames(m_obs )), drop = FALSE]
  pred_high <- m_pred[, grepl(keep_pat, colnames(m_pred)), drop = FALSE]
  if (ncol(obs_high) == 0 || ncol(pred_high) == 0) {
    stop("No columns matched '^(hit|cr)-k[123]-high[12]$'.")
  }
  
  # Extract condition names
  make_denoms <- function(nms) {
    k_num  <- as.integer(sub(".*-k([123])-.*", "\\1", nms))
    is_hit <- grepl("^hit", nms)
    den_hit <- c(65, 50, 35) # k1,k2,k3 targets
    den_cr  <- c(35, 50, 65) # k1,k2,k3 lures
    ifelse(is_hit, den_hit[k_num], den_cr[k_num])
  }
  den_obs  <- make_denoms(colnames(obs_high))
  den_pred <- make_denoms(colnames(pred_high))
  
  # Proportions 
  obs_prop  <- sweep(obs_high,  2, den_obs,  "/")
  pred_prop <- sweep(pred_high, 2, den_pred, "/")
  
  # Long format 
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
      Response = factor(Response, levels = c("Hit","Correct Rejection")),
      Speed    = factor(Speed, levels = c("Fast","Slow")),
      Source   = factor(Source, levels = c("Observed","Predicted"))
    )

  
  p_both <- ggplot(summary_both, aes(x = K, y = mean_prop)) +
    
    # --- OBSERVED: dots + error ---
    geom_errorbar(
      data = subset(summary_both, Source == "Observed"),
      aes(ymin = ci_lower, ymax = ci_upper, group = Speed),
      width = 0.10,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      data = subset(summary_both, Source == "Observed"),
      aes(shape = Speed),
      size = 3,
      stroke = 1,
      color = "black"
    ) +
    
    # --- PREDICTED: lines ---
    geom_line(
      data = subset(summary_both, Source == "Predicted"),
      aes(linetype = Speed, group = Speed),
      linewidth = 0.6,
      color = "black"
    ) +
    
    facet_wrap(~ Response, nrow = 1) +
    scale_y_continuous(name = "Mean Proportion", limits = y_limits) +
    scale_x_discrete(name = "Target Frequency Condition") +
    
    # FAST = filled, SLOW = open
    scale_shape_manual(values = c("Fast" = 16, "Slow" = 1)) +
    
    # FAST = solid line, SLOW = dashed line
    scale_linetype_manual(values = c("Fast" = "solid", "Slow" = "dashed")) +
    
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
      shape = "Observed (RT)",
      linetype = "Predicted (RT)"
    )
  
  list(summary = summary_both, plot = p_both)
}

res_aj <- plot_obs_pred_high(fit.aj.3clk)
print(res_aj$plot)
ggsave(file.path( fig_dir, "observed_vs_predicted_AJ.png"), res_aj$plot, width = 8, height = 3, dpi = 300)

res_2ht <- plot_obs_pred_high(fit.2ht)
print(res_2ht$plot)
ggsave(file.path( fig_dir, "observed_vs_predicted_2HT.png"), res_2ht$plot, width = 8, height = 3, dpi = 300)

res_sdt <- plot_obs_pred_high(fit.sdt)
print(res_sdt$plot)
ggsave(file.path( fig_dir, "observed_vs_predicted_SDT.png"), res_sdt$plot, width = 8,height = 3, dpi = 300)
