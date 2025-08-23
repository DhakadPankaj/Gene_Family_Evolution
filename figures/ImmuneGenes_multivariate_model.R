library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggrepel)


##########################

setwd("Z:/Immune_gene_evolution/")
df2 <- read.table("dNdS_combined.tsv", sep = "\t", header = TRUE)

groups_to_remove <- unique(df2$group[
  df2$Type == "Immune" & is.na(df2$dNdS)
])

#names(df2)[names(df2) == "Other_class"] <- "Cellular response"
#names(df2)[names(df2) == "Other_pathways"] <- "Cellular defense"
#names(df2)[names(df2) == "JAK.STAT"] <- "JAK/STAT"
#names(df2)[names(df2) == "cGas.Sting"] <- "cGAS-STING"

write.table(df2, file = "summary.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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


df_tr <- df_tr %>% select(all_of(c("HOG", "dNdS", "lambda", "proportion", "Symbol", "Type", "Class", "Path"))) %>%
  unique()

path_colors <- c(
  "Non-immune" = "grey",
  "Toll" = "#2771b7",
  "IMD" = "#c71f26",
  "JAK.STAT" = "#6b3e98",
  "cGAS.Sting" = "#59ab48",
  "RNAi" = "#e8872d",
  "JNK" = "#E7298A",
  "MAPK" = "#e8c61c",
  "Cellular.defense" = "#A6CEE3",
  "Multiple_pathways" = "#666666"
)

custom_theme <- theme(plot.title = element_blank(),
                      plot.margin = margin(20, 20, 20, 20),
                      axis.title = element_text(size=16),
                      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
                      axis.title.x = element_text(size = 14, margin = margin(r = 15)),
                      axis.text = element_text(size=11),
                      panel.background = element_rect(fill='transparent'),
                      plot.background = element_rect(fill='transparent', color=NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black" ))


df_tr$Path <- factor(
  df_tr$Path,
  levels = c("Toll", "IMD", "RNAi", "cGAS.Sting", "JAK.STAT", "JNK", "MAPK", "Cellular.defense", "Multiple", "Non-immune"))

#remove raws with NA in dNdS or lambda
df_tr <- df_tr[complete.cases(df_tr[, c("dNdS", "lambda")]), ]

dNdS_threshold <- quantile(subset(df_tr, Type == "Immune")$dNdS, 0.975, na.rm = TRUE)
prop_threshold <- quantile(subset(df_tr, Type == "Immune")$proportion, 0.975, na.rm = TRUE)
lambda_threshold <- quantile(subset(df_tr, Type == "Immune")$lambda, 0.975, na.rm = TRUE)

df_top2.5 <- df_tr %>%
  mutate(
    top_2.5_dNdS = ifelse(Type == "Immune" & dNdS >= dNdS_threshold, "Yes", "No"),
    top_2.5_lambda = ifelse(Type == "Immune" & lambda >= lambda_threshold, "Yes", "No"),
    top_2.5_proportions = ifelse(Type == "Immune" & proportion >= prop_threshold, "Yes", "No")
  )

top_genes <- df_top2.5 %>%
  filter(top_2.5_dNdS == "Yes" | top_2.5_lambda == "Yes" | top_2.5_proportions == "Yes")

write.table(top_genes, file = "top_2.5_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

p1 <- ggplot(df_tr, aes(x = proportion, y= dNdS)) +
  geom_point(aes(colour = Path, alpha = Type), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray", alpha = 0.3) +
  geom_hline(yintercept = dNdS_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_vline(xintercept = prop_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_text_repel(data = subset(df_tr, Type == "Immune"), 
                  aes(label = Symbol), 
                  size = 3,
                  max.overlaps = 10,
                  box.padding = 0.3) +
  scale_color_manual(values = path_colors) +
  scale_alpha_manual(values = c("Non-immune" = 0.2, "Immune" = 0.8)) +
  guides(alpha = "none") +
  custom_theme +
  theme(legend.position = "none")

p1

p2 <- ggplot(df_tr, aes(x = lambda, y= dNdS)) +
  geom_point(aes(colour = Path, alpha = Type), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray", alpha = 0.3) +
  geom_hline(yintercept = dNdS_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_vline(xintercept = lambda_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_text_repel(data = subset(df_tr, Type == "Immune"), 
                  aes(label = Symbol), 
                  size = 3,
                  max.overlaps = 10,
                  box.padding = 0.3) +
  scale_color_manual(values = path_colors) +
  scale_alpha_manual(values = c("Non-immune" = 0.2, "Immune" = 0.8)) +
  guides(alpha = "none") +
  custom_theme

p2
ggsave("legends.svg", get_legend(p2))


p3 <- ggplot(df_tr, aes(x = proportion, y= lambda)) +
  geom_point(aes(colour = Path, alpha = Type), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray", alpha = 0.3) +
  geom_hline(yintercept = lambda_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_vline(xintercept = prop_threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  geom_text_repel(data = subset(df_tr, Type == "Immune"), 
                  aes(label = Symbol), 
                  size = 3,
                  max.overlaps = 10,
                  box.padding = 0.3) +
  scale_color_manual(values = path_colors) +
  scale_alpha_manual(values = c("Non-immune" = 0.2, "Immune" = 0.8)) +
  guides(alpha = "none") +
  custom_theme +
  theme(legend.position = "none")

ggsave("Rates_correlation.svg", plot = plot_grid(p1, p2 + theme(legend.position = "none"), p3, ncol = 2), width = 8.27, height = 11.69/1.5, units = "in")



#######################################################
############# Violin/box plots ########################
#######################################################

### log transform ###
df_tr <- df2 %>%
  select(-c(Receptor, Signalling, Effector, Antiviral, Unclassified, 
            Toll, IMD, JAK.STAT, cGAS.Sting, RNAi, JNK, MAPK, Cellular.defense)) %>%
  left_join(immune_class, by = "HOG") %>%
  mutate(Class = ifelse(Type != "Immune", "Non-immune", Class)) %>%
  left_join(immune_path, by = "HOG") %>%
  mutate(Path = ifelse(Type != "Immune", "Non-immune", Path))

df_tr$log.dNdS<-log(df_tr$dNdS)
df_tr$log.prop <- log(df_tr$proportion)
df_tr$log.prop[df_tr$log.prop<(-8)]<-NA
df_tr$busted_positive <- factor(df_tr$pvalue < 0.001, levels = c(FALSE, TRUE), labels = c("no", "yes"))

hist(df_tr$log.prop)
hist(df_tr$log.dNdS)

df_tr$log.lambda=log(df_tr$lambda)
hist(df_tr$log.lambda)
df_tr$log.lambda[df_tr$log.lambda<(-15)]<-NA
df_tr$is.variable<-vector(length=nrow(df_tr))
df_tr$is.variable[!is.na(df_tr$log.lambda)]<-TRUE

levels(df_tr$Class)
levels(df_tr$Path)

df_plot <- df_tr %>%
  select(c(HOG, Type, Name, log.dNdS, lambda, log.lambda, Class, Path, log.prop, pvalue, busted_positive, group)) %>%
  unique()



#df_plot$log.lambda <- log(df_plot$lambda)

###### Raw dNdS and lambda plots #######
########################################
custom_theme <- theme(plot.title = element_blank(),
                      plot.margin = margin(20, 20, 20, 20),
                      axis.title = element_text(size=16),
                      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
                      axis.title.x = element_blank(),
                      axis.text = element_text(size=11),
                      panel.background = element_rect(fill='transparent'),
                      plot.background = element_rect(fill='transparent', color=NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black" ))

class_colors <- c(
  "Non-immune" = "grey",
  "Antiviral" = "#E41A1C",
  "Effector" = "#377EB8",
  "Unclassified" = "#4DAF4A",
  "Receptor" = "#984EA3",
  "Signalling" = "#FF7F00",
  "Multiple" = "#A65628"
)

path_colors <- c(
  "Non-immune" = "grey",
  "Toll" = "#2771b7",
  "IMD" = "#c71f26",
  "JAK.STAT" = "#6b3e98",
  "cGAS.Sting" = "#59ab48",
  "RNAi" = "#e8872d",
  "JNK" = "#E7298A",
  "MAPK" = "#e8c61c",
  "Cellular.defense" = "#A6CEE3",
  "Multiple_pathways" = "#666666"
)


df_plot$Type <- factor(
  df_plot$Type,
  levels = c("Immune","Non-immune"))

df_plot$Class <- factor(
  df_plot$Class,
  levels = c("Receptor", "Signalling", "Effector", "Antiviral", "Unclassified",  "Multiple", "Non-immune"))

df_plot$Path <- factor(
  df_plot$Path,
  levels = c("Toll", "IMD", "RNAi", "cGAS.Sting", "JAK.STAT", "JNK", "MAPK", "Cellular.defense", "Multiple", "Non-immune"))

#df_plot <- df_plot %>%
#  mutate(log.lambda = ifelse(log.lambda < -15, -15, log.lambda))

df_plot <- df_plot %>%
  mutate(log.dNdS = ifelse(log.dNdS < -5, -5, log.dNdS))

df_plot <- df_plot %>%
  mutate(log.prop = ifelse(log.prop < -6, -6, log.prop))

#metric <- "log.lambda"

create_significance_plot <- function(p_val, plot){
  sig_comparisons <- p_val %>%
    filter(level1 == "Non-immune" & p_adj < 0.05) %>%
    mutate(
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**", 
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  if (nrow(sig_comparisons) > 0) {
    
    # Calculate y-position for asterisks (above the highest point)
    y_max <- max(df_p[[metric]], na.rm = TRUE)
    y_range <- diff(range(df_p[[metric]], na.rm = TRUE))
    asterisk_y <- y_max + 0.05 * y_range
    
    # Create annotation data for asterisks
    asterisk_data <- sig_comparisons %>%
      select(level2, significance) %>%
      rename(Class = level2) %>%
      mutate(y_pos = asterisk_y)
    
    # Add asterisks to the plot
    plot <- plot +
      geom_text(
        data = asterisk_data,
        aes(x = Class, y = y_pos, label = significance),
        color = "black",
        size = 6,
        vjust = 0.5,
        inherit.aes = FALSE
      )
  }
  
  return(plot)
  
}




for (metric in c("log.dNdS", "log.lambda", "log.prop")) {
  plot_title <- if (metric == "log.dNdS") {
    "dN/dS ratio"
  } else if (metric == "log.lambda") {
    "Rate of gene turnover"
  } else if (metric == "log.prop") {
    "Sites under diversifying selection"
  }
  

groups_to_remove <- unique(df_plot$group[
    df_plot$Type == "Immune" & is.na(df_plot$metric)
  ])
#Remove all rows with groups (those contain "NA")
df_p <- df_plot[!df_plot$group %in% groups_to_remove, ]
df_p <- df_p[complete.cases(df_p[,c(metric)]), ]

df_p <- df_p %>% select(all_of(c("HOG", metric, "Type", "Class", "Path"))) %>%
  group_by(HOG, Type) %>%
  summarise(
    across(all_of(metric), ~mean(.x, na.rm = TRUE)),   # mean for metric
    across(c(Class, Path), ~first(.x)),                # keep first value for Class and Path
    n = n(),                                           # count of rows per group
    .groups = "drop"
  )


yvals <- df_p[[metric]]
breaks_raw <- pretty(range(yvals, na.rm = TRUE))
breaks_exp <- unique(round(exp(breaks_raw), 2))
# Remove duplicate labels (e.g., multiple zeros)
breaks_exp <- breaks_exp[!duplicated(breaks_exp)]


if (metric == "log.lambda") {
  breaks_exp <- c(0.001, 0.005, 0.01, 0.02)
}

plot1 <- ggplot(df_p, aes(x = Type, y = .data[[metric]], color = Type)) +
  geom_violin( fill = NA,width = 0.5,  size = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.02) +
  ylab(plot_title) +
  scale_color_manual(values = c("red", "grey")) +
  scale_y_continuous(
    breaks = log(breaks_exp),  
    labels = breaks_exp        
  ) +
  custom_theme +
  theme(legend.position = "none",
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1))

plot1

assign(paste0("plot1", "_", metric), plot1)

p_val <- read.table(paste0("models/", metric, "_pairwise_Class.tsv"), sep= "\t", header=TRUE)

plot2 <- ggplot(df_p, aes(x = Class, y = .data[[metric]], color = Class)) +
  geom_violin( fill = NA, size = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.1) +         
  ylab(plot_title) +
  scale_color_manual(values = class_colors) +
  scale_y_continuous(
    breaks = log(breaks_exp),  
    labels = breaks_exp        
  ) +
  custom_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1))

#plot2

plot2 <- create_significance_plot(p_val, plot2)

assign(paste0("plot2", "_", metric), plot2)

p_val <- read.table(paste0("models/", metric, "_pairwise_Path.tsv"), sep= "\t", header=TRUE)

plot3 <- ggplot(df_p, aes(x = Path, y = .data[[metric]], color = Path)) +
  geom_violin( fill = NA,  size = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.1) +         
  ylab(plot_title) +
  scale_color_manual(values = path_colors) +
  scale_y_continuous(
    breaks = log(breaks_exp),  
    labels = breaks_exp        
  ) +
  custom_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1))

#plot3

plot3 <- create_significance_plot(p_val, plot3)

assign(paste0("plot3", "_", metric), plot3)

}

ggsave("Immune_Vs_Control.svg", plot = plot_grid(plot1_log.dNdS, plot1_log.prop, plot1_log.lambda, ncol = 3), width = 8.27, height = 11.69/3.5, units = "in")

ggsave("Immune_Classes.svg", plot = plot_grid(plot2_log.dNdS, plot2_log.prop, plot2_log.lambda, ncol = 1), width = 8.27/1.2, height = 11.69/1.1, units = "in")

ggsave("Immune_Pathways.svg", plot = plot_grid(plot3_log.dNdS,plot3_log.prop, plot3_log.lambda, ncol = 1), width = 8.27/1.2, height = 11.69/1.1, units = "in")
#######################################################