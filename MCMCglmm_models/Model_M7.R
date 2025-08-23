library(MCMCglmm)
library(ape)
library(dplyr)
library(purrr)
library(tidyr)

#setwd("Z:/Immune_gene_evolution/")

df <- read.table("lambda_dNdS_combined.tsv", sep = "\t", header = TRUE)

#df$Type[df$Type == "Control"] <- "Non-immune"
#names(df)[names(df) == "Other_class"] <- "Unclassified"
#names(df)[names(df) == "Other_pathways"] <- "Cellular defense"
#names(df)[names(df) == "JAK.STAT"] <- "JAK/STAT"
#names(df)[names(df) == "cGas.Sting"] <- "cGAS-STING"

df <- df[complete.cases(df[, c("dNdS", "lambda", "Type", "proportion" ,"median_length", "FPKM", "RSA", "inter", "group")]), ]


df$log.lambda=log(df$lambda)
#hist(df$log.lambda)
df$log.lambda[df$log.lambda<(-15)]<-NA
df$is.variable<-vector(length=nrow(df))
df$is.variable[!is.na(df$log.lambda)]<-TRUE
df$log.dNdS<-log(df$dNdS+0.0001)

df$log.prop <- log(df$proportion)
df$log.prop[df$log.prop<(-8)]<-NA

# FDR correction
#df <- df %>%
#  mutate(
#    p_adj_BH = p.adjust(pvalue, method = "BH"),
#  )

#sum(df$p_adj_ < 0.1)
df$busted_positive <- ifelse(df$pvalue < 0.001, TRUE, FALSE)

#hist(df$log.prop)
#hist(df$log.dNdS)

df <- df %>%
  arrange(FPKM) %>%
  mutate(FPKM_category = ntile(FPKM, 40))

df$Type <- relevel(factor(df$Type), ref = "Non-immune")

#Hadfield approved
prior <- list(
  B = list(mu = rep(0, 30), V = diag(30) * 1e+10),  # Prior for fixed effects
  G = list(G1 = list(V = diag(5), nu = 5, alpha.mu=rep(0,5), alpha.V=diag(5)*1000)),  #Parameter expanded
  R = list(V = diag(5), fix=5, nu = 2.002)  # Prior for residual
)

model <- MCMCglmm(
  cbind(log.dNdS, log.lambda, is.variable, log.prop, busted_positive) ~ trait - 1 + trait:Type + trait:median_length + trait:FPKM_category + trait:RSA + trait:inter,
  random = ~us(trait):group,
  rcov = ~us(trait):units,
  family = c("gaussian", "gaussian","threshold", "gaussian","threshold"),
  prior = prior,
  data = df,
  nitt = 2100000,
  burnin = 100000,
  thin = 100,
  pr = TRUE
)

summ <- summary(model)
sink("multivariate_model_summary.txt")
print(summ)
sink()

#summary(model)

#windows(width = 10, height = 7)
#plot(model)
#par(ask = FALSE)


#####
predictors <- c("Type", "median_length", "FPKM_category", "RSA", "inter")
log_traits <- c("log.dNdS", "log.lambda", "log.prop")
probit_traits <- c("is.variable", "busted_positive")

results_list <- list()

# Function for log traits (proportional change)
process_log_trait <- function(trait, predictor, df) {
    if (predictor == "Type") {
    level <- "Immune"      
    results <- list()
      name <- paste0("trait", trait, ":", predictor, level)
      if (name %in% colnames(model$Sol)) {
        beta <- model$Sol[, name]
        prop_change <- exp(beta) - 1
        
        return(data.frame(
          trait = trait,
          predictor = paste(predictor, level, sep = "_"),
          baseline_prob = NA,
          mean_change = mean(prop_change),
          lower = HPDinterval(as.mcmc(prop_change))[1],
          upper = HPDinterval(as.mcmc(prop_change))[2],
          pMCMC = 2 * min(mean(prop_change > 0), mean(prop_change < 0)),
          interpretation = "proportional_change",
          stringsAsFactors = FALSE
        ))
      }
  } else {
  name <- paste0("trait", trait, ":", predictor)
  
  if (name %in% colnames(model$Sol)) {
    beta <- model$Sol[, name]
    scale_factor <- ifelse(predictor == "median_length", 1000, ifelse(predictor == "RSA", sd(df$RSA, na.rm = TRUE), 1))
    scaled_beta <- beta * scale_factor
    prop_change <- exp(scaled_beta) - 1
    
    return(data.frame(
      trait = trait,
      predictor = predictor,
      baseline_prob = NA,
      mean_change = mean(prop_change),
      lower = HPDinterval(as.mcmc(prop_change))[1],
      upper = HPDinterval(as.mcmc(prop_change))[2],
      pMCMC = 2 * min(mean(prop_change > 0), mean(prop_change < 0)),
      interpretation = "proportional_change"
    ))
  }
  return(NULL)
}

}

# Function for probit traits (probability change)
process_probit_trait <- function(trait, predictor, df) {
    # Handle categorical variables
    if (predictor == "Type") {
      level <- c("Immune")
    
    intercept_name <- paste0("trait", trait)
    
    if (intercept_name %in% colnames(model$Sol)) {
      intercept <- model$Sol[, intercept_name]
        beta_name <- paste0("trait", trait, ":", predictor, level)
        if (beta_name %in% colnames(model$Sol)) {
          beta <- model$Sol[, beta_name]
          
          # Baseline probability (reference category)
          p0 <- pnorm(intercept)
          # Probability for this category
          p1 <- pnorm(intercept + beta)
          prob_change <- p1 - p0
          
          return(data.frame(
            trait = trait,
            predictor = paste(predictor, level, sep = "_"),
            baseline_prob = mean(p0),
            mean_change = mean(prob_change),
            lower = HPDinterval(as.mcmc(prob_change))[1],
            upper = HPDinterval(as.mcmc(prob_change))[2],
            pMCMC = 2 * min(mean(prob_change > 0), mean(prob_change < 0)),
            interpretation = "probability_change",
            stringsAsFactors = FALSE
          ))
        }
    }
  } else {
  beta_name <- paste0("trait", trait, ":", predictor)
  intercept_name <- paste0("trait", trait)
  
  if (beta_name %in% colnames(model$Sol) && intercept_name %in% colnames(model$Sol)) {
    beta <- model$Sol[, beta_name]
    intercept <- model$Sol[, intercept_name]
    
    # Scale factor for meaningful change
    scale_factor <- ifelse(predictor == "median_length", 1000, ifelse(predictor == "RSA", sd(df$RSA, na.rm = TRUE), 1))

    # Baseline probability
    p0 <- pnorm(intercept)
    
    # Probability after change
    p1 <- pnorm(intercept + beta * scale_factor)
    
    # Change in probability
    prob_change <- p1 - p0
    
    return(data.frame(
      trait = trait,
      predictor = predictor,
      baseline_prob = mean(p0),
      mean_change = mean(prob_change),
      lower = HPDinterval(as.mcmc(prob_change))[1],
      upper = HPDinterval(as.mcmc(prob_change))[2],
      pMCMC = 2 * min(mean(prob_change > 0), mean(prob_change < 0)),
      interpretation = "probability_change"
    ))
  }
  return(NULL)
}

}

# Process log traits
for (trait in log_traits) {
  for (pred in predictors) {
    result <- process_log_trait(trait, pred, df)
    if (!is.null(result)) {
      results_list[[paste(trait, pred, sep = "_")]] <- result
    }
  }
}

# Process probit traits
for (trait in probit_traits) {
  for (pred in predictors) {
    result <- process_probit_trait(trait, pred, df)
    if (!is.null(result)) {
      results_list[[paste(trait, pred, sep = "_")]] <- result
    }
  }
}

# Combine results
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

write.table(results_df, "mcmcglmm_results.tsv",sep = "\t" ,row.names = FALSE, quote = FALSE)
