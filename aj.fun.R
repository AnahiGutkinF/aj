#####===========================================================================
##### [crrm] confidence and rt recognition models
##### Script for functions
##### Author: Anahi Gutkin
##### Date: 23/09/2022
##### Version: v5 (04/12/2024)
#####===========================================================================
#-------------------------------------------------------------------------------
# 1. Instal packages and load them.
#-------------------------------------------------------------------------------

if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", dependencies = TRUE)
  library("data.table")
}

#-------------------------------------------------------------------------------
# 2. Function for discretizing continuous variables (y) in bins
#-------------------------------------------------------------------------------

recateg.2bin <- function(data, NSets=47, criteria="geom", groupby = "x"){
  S <- NSets
  L <- length(levels(data[, groupby]))
  frequencies <- matrix(NA, nrow = S, ncol = L*2 + 1, byrow = F)
  for(s in 1:S){
    recog <- data[data$id==s,]
    if (criteria=="median"){
      y_mean <- median(recog$y)
    } else {
      if (criteria=="mean"){
        y_mean <- mean(recog$y)
      } else {
        y_mean <- exp(mean(log(recog$y)))
      }
    }
    recog$y_cat <- as.numeric(recog$y > y_mean)
    mptrt_categories <- paste0( rep(levels(data[, groupby]), each=2), c("0", "1"))
    recog$mptrt_category <- paste0(recog[, groupby],recog$y_cat)
    recog$mptrt_category <- factor(recog$mptrt_category,levels =  mptrt_categories)
    frequencies[s,]<- c( table(recog$mptrt_category), y_mean)
  }
  colnames(frequencies) <- c(mptrt_categories, "break")
  return(frequencies)
}

#-------------------------------------------------------------------------------
# 3. Change factor levels names
#-------------------------------------------------------------------------------

change_factor_levels <- function(factor_var, new_levels) {
  factor_var <- factor(factor_var)
  levels(factor_var) <- new_levels
  return(factor_var)
}

#-------------------------------------------------------------------------------
# 4. Compare nested models
#-------------------------------------------------------------------------------

test.nested <- function(model.gen, model.rest) {
  G2 = round(-2 * (model.rest$goodness.of.fit$individual$Log.Likelihood - 
                     model.gen$goodness.of.fit$individual$Log.Likelihood), 2)
  G2 = ifelse(G2 < 0, 0, G2)
  df = model.gen$model.info$individual$n.parameters - model.rest$model.info$individual$n.parameters
  p.val = round(pchisq(G2, df = 1, lower.tail = FALSE), 3)
  return(data.frame(G2, p.val, df))
}

#-------------------------------------------------------------------------------
# 5. Calculate AIC weights
#-------------------------------------------------------------------------------

calculate_wAIC <- function(muestras,
                           model.names) {
  AIC_values <- muestras[[1]]$information.criteria$individual$AIC
  for (i in 2:length(muestras)) {
    AIC2 <- muestras[[i]]$information.criteria$individual$AIC
    AIC_values <- cbind.data.frame(AIC_values, AIC2)
  }
  colnames(AIC_values) <- model.names
  waic <- function(AIC_val){delta_AIC <- AIC_val - min(AIC_val)
  exponents <- exp(-0.5 * delta_AIC)
  wAIC_weights <- exponents / sum(exponents)
  return(round(wAIC_weights,3))
  }
  wAIC <- as.data.frame(t(apply(AIC_values, 1, FUN = waic)))
  return(wAIC)
}

#-------------------------------------------------------------------------------
# 6. Change from full continuous format to freq
#-------------------------------------------------------------------------------

change_to_freq <- function(datos){
  Cond <- max(datos$Target_frequency)
  I <- max(datos$Id)
  # Cambiar los niveles de los factores
  datos$Stimulus <- change_factor_levels(as.factor(datos$Stimulus), c("learned_item","new_item" ))
  datos$Type_of_trial <- change_factor_levels(as.factor(datos$Type_of_trial), c("hit", "miss", "fa", "cr"))
  datos$Confidence_level <- change_factor_levels(as.factor(datos$Confidence_level), c("high", "low"))
  
  l.x <- paste0(rep(c("hit", "miss", "fa", "cr"), Cond), "-k", rep(1:Cond, each = 4))
  l.cnf <- paste0(rep(l.x, each = 2), c("-high", "-low"))
  datos <- data.frame(
    tree   = datos$Stimulus,
    x      = factor(paste0(datos$Type_of_trial, "-k", datos$Target_frequency), levels = l.x),
    x_cnf  = factor(paste0(datos$Type_of_trial, "-k", datos$Target_frequency, "-", datos$Confidence_level), levels = l.cnf),
    y      = datos$RT,
    id     = datos$Id
  )
  
  # Recategorización y combinación de resultados
  d.cnf.rt <- cbind(recateg.2bin(data = datos, NSets = I, criteria = "geom",
                                 groupby = "x_cnf"))
  return(d.cnf.rt[,1:(ncol(d.cnf.rt)-1)])
}

