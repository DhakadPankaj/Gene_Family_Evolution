library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(purrr)
library(readr)

setwd("Z:/Immune_gene_evolution/HyPhy_MEME_FEL/")

info <- read.table("HOG_RNA_info.tsv", sep = "\t", header = TRUE)
info2 <- read_excel("MEME_FEL_results.xlsx", sheet = 1)
info2 <- info2[,c("HOG","Class")]
info <- merge(info, info2, by = "HOG")

classes <- unique(info$Class)

# Domains data (clean it)
domains <- read_excel("conserved_domains.xlsx", sheet = 1)
domains$Query <- str_extract(domains$Query, "rna-[^\\s]+")

domains <-  domains[complete.cases(domains[, "Query"]), ]

# MEME Sites
data_list <- list.files("sites_mappings_MEME", pattern = "*.tsv", full.names = TRUE) %>%
  set_names(basename(.)) %>%
  map(~read_tsv(.x))

# Function to create gene plot
plot_gene_class <- function(class_name, info_data, domains_data, meme_data_list) {
  
  # Get genes for this class
  class_genes <- info_data %>% filter(Class == class_name)
  
  if(nrow(class_genes) == 0) return(NULL)
  
  # Pre-filter genes that have MEME sites
  genes_with_meme <- c()
  
  for(i in 1:nrow(class_genes)) {
    gene_info <- class_genes[i, ]
    hog_id <- gene_info$HOG
    
    # Check if this gene has MEME sites
    meme_file_pattern <- paste0(hog_id, "_site_mappings.tsv")
    matching_files <- names(meme_data_list)[grepl(meme_file_pattern, names(meme_data_list))]
    
    if(length(matching_files) > 0) {
      meme_sites <- filter(meme_data_list[[matching_files[1]]], 
                           Sequence_ID == gene_info$RNA & Unmasked_Position != "N/A")
      
      if(nrow(meme_sites) > 0) {
        positions <- as.numeric(meme_sites$Unmasked_Position)
        if(length(positions[!is.na(positions)]) > 0) {
          genes_with_meme <- c(genes_with_meme, i)
        }
      }
    }
  }
  
  # Filter class_genes to only include those with MEME sites
  if(length(genes_with_meme) == 0) {
    cat("No genes with MEME sites found in class:", class_name, "\n")
    return(NULL)
  }
  
  class_genes <- class_genes[genes_with_meme, ]
  cat("Plotting", nrow(class_genes), "genes with MEME sites in class:", class_name, "\n")
  
  domains_data <- domains_data %>%
    filter(Query %in% class_genes$RNA)
  
  unique_domains <- unique(domains_data$`Short name`)
  # Use light colors - you can customize this palette
  domain_colors <- RColorBrewer::brewer.pal(min(length(unique_domains), 8), "Set2")
  if(length(unique_domains) > 8) {
    domain_colors <- rainbow(length(unique_domains), alpha = 0.4)
  }
  names(domain_colors) <- unique_domains
  
  # Prepare plotting data
  plot_data <- data.frame()
  domain_data <- data.frame()
  meme_sites_data <- data.frame()
  
  class_genes <- class_genes[order(class_genes$Symbol), ]
  
  for(i in 1:nrow(class_genes)) {
    gene_info <- class_genes[i, ]
    gene_length <- gene_info$Length  
    gene_name <- gene_info$Symbol    
    hog_id <- gene_info$HOG
    
    # Add gene backbone
    gene_backbone <- data.frame(
      gene = gene_name,
      start = 0,
      end = gene_length,
      y_pos = i,
      type = "gene"
    )
    plot_data <- rbind(plot_data, gene_backbone)
    
    # Add domains for this gene
    if(!is.null(gene_info$RNA) && !is.na(gene_info$RNA)) {
      gene_domains <- domains_data %>% 
        filter(Query == gene_info$RNA) %>% 
        select(From, To, domain_name = `Short name`) %>% 
        mutate(
          gene = gene_name,
          y_pos = i,
          type = "domain"
        )
      
      if(nrow(gene_domains) > 0) {
        names(gene_domains)[1:2] <- c("start", "end")
        domain_data <- rbind(domain_data, gene_domains)
      }
    }
    
    # Add MEME sites
    meme_file_pattern <- paste0(hog_id, "_site_mappings.tsv")
    matching_files <- names(meme_data_list)[grepl(meme_file_pattern, names(meme_data_list))]
    
    if(length(matching_files) > 0) {
      meme_sites <- filter(meme_data_list[[matching_files[1]]], Sequence_ID==gene_info$RNA  & Unmasked_Position != "N/A")
      if("Unmasked_Position" %in% names(meme_sites)) {
        meme_positions <- data.frame(
          gene = gene_name,
          position =  as.numeric(meme_sites$Unmasked_Position),
          y_pos = i,
          type = "meme_site"
        )
        meme_sites_data <- rbind(meme_sites_data, meme_positions)
      }
    }
  }
  
  # Create the plot
  p <- ggplot() +
    # Gene backbone (white/light grey)
    geom_rect(data = plot_data, 
              aes(xmin = start, xmax = end, ymin = y_pos - 0.2, ymax = y_pos + 0.2),
              fill = "lightgrey", color = "black", size = 0.3) +
    
    # Domains (grey)
    {if(nrow(domain_data) > 0) {
      geom_rect(data = domain_data,
                aes(xmin = start, xmax = end, ymin = y_pos - 0.2, ymax = y_pos + 0.2, fill = domain_name),
                color = "black", size = 0.3)
    }} +
    
    {if(nrow(domain_data) > 0) {
      scale_fill_manual(values = domain_colors, name = "Domain")
    }} +
    
    # Domain labels
    {if(nrow(domain_data) > 0) {
      geom_text(data = domain_data %>% 
                  mutate(domain_width = end - start,
                         use_repel = domain_width < 50),  # Repel for short domains
                aes(x = (start + end) / 2, y = y_pos + 0.42, label = domain_name),
                size = 3.5,
                hjust = 0.5,
                vjust = 0.5,
                check_overlap = TRUE)  # Automatically hide overlapping text
    }} +
    
    # MEME sites (red points)
    {if(nrow(meme_sites_data) > 0) {
      geom_point(data = meme_sites_data,
                 aes(x = position, y = y_pos),
                 color = "brown", size = 3.5, shape = "|")
    }} +
    
    # Gene labels
    geom_text(data = plot_data,
              aes(x = -max(end) * 0.01, y = y_pos, label = gene),
              hjust = 1, vjust = 0.5, size = 5) +
    geom_text(data = plot_data,
              aes(x = end + max(end) * 0.01, y = y_pos, label = paste0(end, " AA")),
              hjust = 0, vjust = 0.5, size = 3.5, color = "grey30") +
    # Formatting
    scale_y_continuous(breaks = 1:nrow(class_genes), 
                       labels = class_genes$Symbol,
                       limits = c(0.5, nrow(class_genes) + 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    
    labs(
      title = paste(class_name),
      x = "",
      y = ""
    ) +
    
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )
  return(p)
}
p
# Generate plots for each class
plots <- list()
for (class in classes) {
  cat("Processing class:", class, "\n")
  p <- plot_gene_class(class, info, domains, data_list)
  
  if(!is.null(p)) {
    plots[[class]] <- p
    
    ggsave(paste0("gene_structure_", gsub("[^A-Za-z0-9]", "_", class), ".svg"), plot = p, 
           width = 8.27, height = nrow(info[info$Class == class, ]) * 0.5, units = "in")
    
    # Display plot
    print(p)
  }
}

# plots[["class_name"]]
