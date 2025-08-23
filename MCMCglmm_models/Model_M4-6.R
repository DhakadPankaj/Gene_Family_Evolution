library(MCMCglmm)
library(ape)
library(dplyr)
library(tidyr)
library(parallel)

#library(ggplot2)
#library(cowplot)

##########################

setwd("/home/s2215768/Immune_gene_evolution/models")
df2 <- read.table("dNdS_combined.tsv", sep = "\t", header = TRUE)

groups_to_remove <- unique(df2$group[
  df2$Type == "Immune" & is.na(df2$dNdS)
])

#Remove all rows with groups (those contain "NA")
df2 <- df2[!df2$group %in% groups_to_remove, ]
df2 <- df2[complete.cases(df2[,c("dNdS")]), ]

# For Immune rows
immune_class <- df2 %>%
  filter(Type == "Immune") %>%
  pivot_longer(
    cols = c(Receptor, Signalling, Effector, Antiviral, Unclassified),
    names_to = "Class", values_to = "Presence"
  ) %>%
  filter(Presence == "Yes") %>%
  group_by(HOG) %>%
  summarise(
    Class = if (n() == 1) first(Class) else "Multiple",
    .groups = "drop"
  )

immune_path <- df2 %>%
  filter(Type == "Immune") %>%
  pivot_longer(
    cols = c(Toll, IMD, JAK.STAT, cGAS.Sting, RNAi, JNK, MAPK, Cellular.defense),
    names_to = "Path", values_to = "Presence"
  ) %>%
  filter(Presence == "Yes") %>%
  group_by(HOG) %>%
  summarise(
    Path = if (n() == 1) first(Path) else "Multiple",
    .groups = "drop"
  )

# Combine both
df_tr <- df2 %>%
  select(-c(Receptor, Signalling, Effector, Antiviral, Unclassified, 
            Toll, IMD, JAK.STAT, cGAS.Sting, RNAi, JNK, MAPK, Cellular.defense)) %>%
  left_join(immune_class, by = "HOG") %>%
  mutate(Class = ifelse(Type != "Immune", "Non-immune", Class)) %>%
  left_join(immune_path, by = "HOG") %>%
  mutate(Path = ifelse(Type != "Immune", "Non-immune", Path))


df_tr$Class <- as.factor(df_tr$Class)
df_tr$Path <- as.factor(df_tr$Path)
df_tr$group <- as.factor(df_tr$group)

### log transform ###
df_tr$log.dNdS<-log(df_tr$dNdS+0.0001)
df_tr$log.prop <- log(df_tr$proportion)
df_tr$log.prop[df_tr$log.prop<(-8)]<-NA
df_tr$busted_positive <- factor(df_tr$pvalue < 0.001, levels = c(FALSE, TRUE), labels = c("no", "yes"))

#hist(df_tr$log.prop)
#hist(df_tr$log.dNdS)

df_tr$log.lambda=log(df_tr$lambda)
#hist(df_tr$log.lambda)
df_tr$log.lambda[df_tr$log.lambda<(-15)]<-NA
df_tr$is.variable<-vector(length=nrow(df_tr))
df_tr$is.variable[!is.na(df_tr$log.lambda)]<-TRUE

levels(df_tr$Class)
levels(df_tr$Path)


### Pairwise comparisons
### Pairwise comparisons
pairwise_compare <- function(trait, predictor, model, levels, back_transform = TRUE) {
  # Get the levels of the predictor from data
  predictor_levels <- levels
  if (is.null(predictor_levels)) stop("Could not determine predictor levels. Check data table")
  
  # Get posterior samples for each level (on linear scale)
  posteriors_linear <- list()
  # Reference level (first in levels)
  posteriors_linear[[predictor_levels[1]]] <- model$Sol[, paste0("trait", trait)]
  for (lev in predictor_levels[-1]) {
    colname <- paste0("trait", trait, ":", predictor, lev)
    posteriors_linear[[lev]] <- model$Sol[, paste0("trait", trait)] + model$Sol[, colname]
  }
  
  # Back-transform if needed
  posteriors_data_scale <- list()
  for (lev in predictor_levels) {
    if (back_transform && grepl("^log\\.", trait)) {
      # For log-transformed traits, exponentiate
      posteriors_data_scale[[lev]] <- exp(posteriors_linear[[lev]])
    } else if (back_transform && trait %in% c("is.variable", "busted_positive")) {
      # For probit traits, convert to probabilities
      posteriors_data_scale[[lev]] <- pnorm(posteriors_linear[[lev]])
    } else {
      # No transformation needed
      posteriors_data_scale[[lev]] <- posteriors_linear[[lev]]
    }
  }
  
  # Pairwise comparisons on data scale
  res <- data.frame(
    trait = character(),
    predictor = character(),
    level1 = character(),
    level2 = character(),
    mean_level1 = numeric(),
    mean_level2 = numeric(),
    mean_diff = numeric(),
    fold_change = numeric(),  # Only for log traits
    lower = numeric(),
    upper = numeric(),
    pMCMC = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(predictor_levels)-1)) {
    for (j in (i+1):length(predictor_levels)) {
      l1 <- predictor_levels[i]
      l2 <- predictor_levels[j]
      
      val1 <- posteriors_data_scale[[l1]]
      val2 <- posteriors_data_scale[[l2]]
      diff <- val1 - val2
      
      mean_val1 <- mean(val1)
      mean_val2 <- mean(val2)
      mean_diff <- mean(diff)
      ci <- quantile(diff, c(0.025, 0.975))
      pval <- 2 * min(mean(diff > 0), mean(diff < 0))
      
      # Calculate fold change for log traits
      fold_change <- ifelse(grepl("^log\\.", trait), mean_val1 / mean_val2, NA)
      
      res <- rbind(res, data.frame(
        trait = trait,
        predictor = predictor,
        level1 = l1,
        level2 = l2,
        mean_level1 = mean_val1,
        mean_level2 = mean_val2,
        mean_diff = mean_diff,
        fold_change = fold_change,
        lower = ci[1],
        upper = ci[2],
        pMCMC = pval
      ))
    }
  }

  write.table(res, paste0(trait, "_pairwise_", predictor, ".tsv"), 
              row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
  
}  

process_log_trait <- function(trait, predictor, model, df) {
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
process_probit_trait <- function(trait, predictor, model, df) {
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



df_tr <- df_tr %>%
  arrange(FPKM) %>%
  mutate(FPKM_category = ntile(FPKM, 40))

df_lambda <- df_tr[complete.cases(df_tr[, c("median_length", "FPKM", "RSA", "inter", "group")]), ]
df_lambda$Type <- relevel(factor(df_lambda$Type), ref = "Non-immune")
df_lambda$Class <- relevel(factor(df_lambda$Class), ref = "Non-immune")
df_lambda$Path <- relevel(factor(df_lambda$Path), ref = "Non-immune")

#Prior for the models
prior4 <- list(
  B = list(mu = rep(0, 12), V = diag(12) * 1e+10),  # Prior for fixed effects
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)),  #Parameter expanded
  R = list(V = diag(2), fix=2, nu = 2.002)  # Prior for residual
)

prior5 <- list(
  B = list(mu = rep(0, 22), V = diag(22) * 1e+10),  # Prior for fixed effects
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)),  #Parameter expanded
  R = list(V = diag(2), fix=2, nu = 2.002)  # Prior for residual
)

prior6 <- list(
  B = list(mu = rep(0, 28), V = diag(28) * 1e+10),  # Prior for fixed effects
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)),  #Parameter expanded
  R = list(V = diag(2), fix=2, nu = 2.002)  # Prior for residual
)

classes <- levels(df_lambda$Class)
pathways <- levels(df_lambda$Path)

#############################################################
############## Gene Turnover ################################

# Run models 4, 5, 6 in parallel
run_model <- function(model_num, df_lambda, prior, formula, filename_prefix, predictor_levels, predictor_name) {
    library(MCMCglmm)
    model <- MCMCglmm(
        formula,
        random = ~us(trait):group,
        rcov = ~us(trait):units,
        family = c("gaussian", "threshold"),
        prior = prior,
        data = df_lambda,
        nitt = 2100000,
        burnin = 100000,
        thin = 100,
        pr = TRUE
    )
    summ <- summary(model)
    sink(paste0(filename_prefix, "_summary.txt"))
    print(summ)
    sink()
  if (!is.null(predictor_levels)) {
        pairwise_compare("log.lambda", predictor_name, model, predictor_levels)
        pairwise_compare("is.variable", predictor_name, model, predictor_levels)
    }

      # back-transform log traits

   if (is.null(predictor_name)) {
    predictors <- c("Type", "median_length", "FPKM_category", "RSA", "inter")
  } else {
    predictors <- c("median_length", "FPKM_category", "RSA", "inter")
  }

  log_traits <- c("log.lambda")
  probit_traits <- c("is.variable")

  results_list <- list()

  for (trait in log_traits) {
  for (pred in predictors) {
    result <- process_log_trait(trait, pred, model, df_lambda)
    if (!is.null(result)) {
      results_list[[paste(trait, pred, sep = "_")]] <- result
    }
  }
}

# Process probit traits
for (trait in probit_traits) {
  for (pred in predictors) {
    result <- process_probit_trait(trait, pred, model, df_lambda)
    if (!is.null(result)) {
      results_list[[paste(trait, pred, sep = "_")]] <- result
    }
  }
}

  # Combine results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  write.table(results_df, paste0(filename_prefix, "_summary_V2.txt"),sep = "\t" ,row.names = FALSE, quote = FALSE)

    return(model)
}

# Prepare arguments for each model
args_list <- list(
    list(
        model_num = 4,
        df_lambda = df_lambda,
        prior = prior4,
        formula = cbind(log.lambda, is.variable) ~ trait - 1 + trait:Type + trait:median_length + trait:FPKM_category + trait:RSA + trait:inter,
        filename = "lambda_Type",
        pairwise = FALSE,
        levels = NULL,
        predictor = NULL
    ),
    list(
        model_num = 5,
        df_lambda = df_lambda,
        prior = prior5,
        formula = cbind(log.lambda, is.variable) ~ trait - 1 + trait:Class + trait:median_length + trait:FPKM_category + trait:RSA + trait:inter,
        filename = "lambda_Class",
        pairwise = TRUE,
        levels = classes,
        predictor = "Class"
    ),
    list(
        model_num = 6,
        df_lambda = df_lambda,
        prior = prior6,
        formula = cbind(log.lambda, is.variable) ~ trait - 1 + trait:Path + trait:median_length + trait:FPKM_category + trait:RSA + trait:inter,
        filename = "lambda_Path",
        pairwise = TRUE,
        levels = pathways,
        predictor = "Path"
    )
)


cl2 <- makeCluster(3)
clusterExport(cl2, varlist = c("run_model", "df_lambda", "prior4", "prior5", "prior6", "classes", "pathways", "pairwise_compare", "process_log_trait", "process_probit_trait"))
clusterEvalQ(cl2, {library(MCMCglmm)})

models <- parLapply(cl2, args_list, function(args) {
  run_model(args$model_num, args$df_lambda, args$prior, args$formula, args$filename, args$levels, args$predictor)
})

stopCluster(cl2)