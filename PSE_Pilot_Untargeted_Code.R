###load libraries####
library(ggplot2) 
library(RColorBrewer)
library(tidyverse)
library(rcompanion)
library(writexl)
library(rstatix)
library(tibble)
library(lmtest)
library(gridExtra)
library(grid)
library(dplyr)
library(ggpubr)
library(tidyr)
library(stringr) 
library(broom)
library(Metrics)
library(devtools)
library(dplyr)
library(gridExtra)
library(cowplot)
library(plyr)
library(naniar)
library(ggh4x)
library(EnvStats)
library(readr)
library(ggvenn)
library(ggVennDiagram)
library(readxl)
library(pheatmap)
library(janitor)
library(ggrepel)
library(SciViews) #to get ln function
library(ape)
library(QSUR)

#before importing anything to R, make sure your data was saved from MSDial as text. Otherwise, there will be rows of data that are misaligned with the columns#
#######write custom function for cleaning ms dial outputs#######
ms_dial_clean <- function(df){
  df <- as.data.frame(df)
  df <- df[-c(1:3),]  #remove unwanted rows#
  names(df) <- df[1,] #make first row the column names#
  df <- df[-c(1),] #remove duplicate first row#
  names(df) <- gsub(" ", "_", names(df)) #replace spaces in column names with underscore#
  names(df) <- gsub("/", ".", names(df)) #replace slashes in column names with period#
  names(df) <- gsub("-", "_", names(df)) #replace dashes in column names with underscore#
  names(df) <- gsub("%", "percent", names(df)) #replace percentages with the word percent#
  index <- match("1",names(df)) #create variable for column index for the first instance of "1"#
  colnames(df)[index] <- "Sample.Avg" #set the corresponding column for sample average from 1 to sample.avg#
  colnames(df)[index+1] <- "Sample.Std" #set the corresponding column for sample std from 1 sto sample.std"#
  df$Sample.Avg = as.numeric(df$Sample.Avg)
  df$Sample.Std = as.numeric(df$Sample.Std)
  df <- df %>% as.data.frame(row.names = 1:nrow(.)) #reset row number#
  df <- subset(df, select = -c(Post_curation_result, Fill_percent, Isotope_tracking_parent_ID,
                               Isotope_tracking_weight_number, Spectrum_reference_file_name, 
                               Weighted_dot_product, Reverse_dot_product,Simple_dot_product)) #remove unneeded columns#
  df #so output is a data frame#
}
lod_and_percent_detection <- function(df){ 
  var_name <- colnames(df[,grep("_S_", colnames(df))]) #set variable for participant samples#
  df[var_name][df[var_name]==0] <- NA #make all 0s NAs in participant samples#
  df[var_name][df[var_name]<0] <- NA #make all 0s NAs in participant samples#
  df[,'Detection_Count'] = NA
  df <- df %>% mutate(Detection_Count = rowSums(!is.na(.[grepl('_S_', colnames(.))]))) #change the value in detection count to the number of columns not equal to NA#
  number_samples_manual = sum(grepl(paste0("_S_"),names(df))) #count number of deployed samples in total#
  df[,'Percent_Detection'] = df$Detection_Count*100/number_samples_manual
  df<-df %>%
    dplyr::mutate(across(all_of(var_name), as.numeric))
  index2 <- match('Detection_Count', names(df)) #create variable for column index for the end of sample columns#
  index3 <- match("Total_score",names(df)) #create variable for column index for the start of sample columns#
  start_sample_column = index3 + 1
  end_sample_column = index2 - 1
  df[, "LOD"] <- apply(df[var_name], 1, min, na.rm = TRUE) #get minimum value across sample columns into a new column, ignoring NAs#
  df <- df %>% 
    mutate(LOD = as.numeric(LOD))%>% #turn the non-zero-min column into numeric#
    mutate(LOD = LOD/5) #get LOD and enter it in the non_zero_min column#
  for (var in var_name){
    df[[var]][is.na(df[[var]])] <- df$LOD[is.na(df[[var]])]}
  df
} 
remove_duplicates <- function(df, grouping, factor){
  MERGED_ALL_NO_DUPS_PER_METHOD_SUBSET_LOD_Corrected <- df %>%
    group_by(SMILES, {{grouping}}) %>% 
    filter(!is.na(SMILES) & SMILES != 'null')%>% #remove anything not matched to a library
    slice_max(order_by = {{factor}}, n =1) %>%#take highest of each unique SMILE and type combination for QC Pool blank sub
    ungroup()
  
  MERGED_ALL_NO_DUPS_PER_METHOD_SUBSET_LOD_Corrected <- MERGED_ALL_NO_DUPS_PER_METHOD_SUBSET_LOD_Corrected %>%
    group_by(INCHIKEY, {{grouping}}) %>% 
    filter(!is.na(INCHIKEY) & INCHIKEY != 'null')%>% 
    slice_max(order_by = {{factor}}, n =1) %>%#take highest of each unique SMILE and type combination for QC Pool blank sub
    ungroup()
  
  df <- MERGED_ALL_NO_DUPS_PER_METHOD_SUBSET_LOD_Corrected %>%
    group_by(Metabolite_name, {{grouping}}) %>% 
    filter(!is.na(Metabolite_name) & Metabolite_name != 'null' & Metabolite_name != 'Unknown')%>% 
    slice_max(order_by = {{factor}}, n =1) %>%#take highest of each unique SMILE and type combination for QC Pool blank sub
    ungroup()
}
blank_correction <- function(df){
  #make separate data frames for blanks and samples to be merged
  PSE_blanks <- df %>% filter(Sample_Type == 'Blank') %>% 
    subset(select = c(Sample, ISTD, Alignment_ID, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, SMILES, m.z_similarity, MS2_matched, Total_score, ISTD_Corrected_Value, Sample_Type, Coating, Location))
  PSE_Samples <- df %>% filter(Sample_Type == 'Sample') %>%
    subset(select = c(Sample, ISTD, Alignment_ID, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, SMILES, m.z_similarity, MS2_matched, Total_score, ISTD_Corrected_Value, Sample_Type, Coating, Location))
  
  #merge blanks with samples, matching by feature, location and material
  
  PSE_Blanks_Matched <- merge(PSE_blanks, PSE_Samples, by = c('Location', 'Coating', 'Alignment_ID', 'Average_Mz', 'Metabolite_name', 'MS.MS_assigned', 'INCHIKEY', 'SMILES', 'm.z_similarity', 'MS2_matched', 'Total_score'), all = T) %>%
    dplyr::rename(Blank_Area = ISTD_Corrected_Value.x) %>%
    dplyr::rename(Sample_Area = ISTD_Corrected_Value.y) %>%
    subset(select = -c(Sample.x, Location, Coating, ISTD.x, ISTD.y, Sample_Type.y, Sample_Type.x))
  
  nrow(PSE_Blanks_Matched) - nrow(PSE_Samples) #sanity check - should give 0
  
  PSE_Blanks_Matched$Sample_Area <- as.numeric(PSE_Blanks_Matched$Sample_Area)
  PSE_Blanks_Matched$Blank_Area <- as.numeric(PSE_Blanks_Matched$Blank_Area)
  
  PSE_Blanks_Matched[,'Blank_Corrected_Area'] = PSE_Blanks_Matched$Sample_Area - PSE_Blanks_Matched$Blank_Area #make blank corrected values a column
  PSE_Blanks_Matched[,'Blank_Sample_Ratio'] = PSE_Blanks_Matched$Blank_Area/PSE_Blanks_Matched$Sample_Area #blank:sample ratio as another column
  
  #pivot wide again
  df <- PSE_Blanks_Matched %>% 
    filter(Blank_Sample_Ratio < 0.35) %>% #filter out features/samples with high blanks
    subset(select = -c(Blank_Area, Sample_Area, Blank_Sample_Ratio))%>% #remove interfering columns
    pivot_wider(names_from = Sample.y, values_from = Blank_Corrected_Area) 
  
  df}
redo_analysis <- function(df){
  df[,'Blank_Sample_Ratio'] = df$Blank_Area/df$Sample_Area
  df[,'Blank_Corrected_Area'] = df$Sample_Area-df$Blank_Area
  
  PSE_Wide_Blank_Sub <- df %>% 
    filter(Blank_Sample_Ratio < 0.35) %>%
    subset(select = -c(Blank_Area, Sample_Area, Blank_Sample_Ratio))%>%
    pivot_wider(names_from = Sample, values_from = Blank_Corrected_Area) 
  
  #######replace with LOD####
  PSE_Wide_Blank_Sub_LODs <- PSE_Wide_Blank_Sub %>% 
    lod_and_percent_detection() %>%
    filter(Percent_Detection > 50 | Percent_Detection == 50) 
  
  PSE_Wide_Blank_Sub_LODs$Sample_Avg =apply(PSE_Wide_Blank_Sub_LODs[grep("JO_S", colnames(PSE_Wide_Blank_Sub_LODs), value=TRUE)],1,mean) #get average across samples
  
  ######remove duplicates - by mz similarity#####
  PSE_Wide_Blank_Sub_LODs$m.z_similarity <- as.numeric(PSE_Wide_Blank_Sub_LODs$m.z_similarity)
  PSE_dup_int1 <- remove_duplicates(df = PSE_Wide_Blank_Sub_LODs, factor = m.z_similarity)
  #####remove duplicates - by sample average#####
  PSE_Dups_removed <- remove_duplicates(df = PSE_dup_int1, factor = Sample_Avg)
  
  ####filter for mz 0.95#####
  PSE_Dups_removed$m.z_similarity <- as.numeric(PSE_Dups_removed$m.z_similarity)
  PSE_mz1_initial <-PSE_Dups_removed %>% filter(m.z_similarity> 0.95) ###NBD note- NA's introduced by coercion warning?
  PSE_remaining_dups <- PSE_mz1_initial %>% group_by(Metabolite_name) %>% filter(n() > 1)
  #view(PSE_remaining_dups) #these all have slightly different average_mz-> pick the higher alignment ID (this is kind of arbitrary)#
  
  
  PSE_remaining_dups_fixed <- PSE_remaining_dups %>%
    group_by(Metabolite_name) %>%
    arrange(Metabolite_name, desc(Alignment_ID)) %>%  # First order by m.z_similarity, then QC_Pooled within each m.z_similarity
    slice_max(order_by = Alignment_ID, n = 1) %>%  # Keep the row with the largest QC_Pooled value for each metabolite group's largest mz
    ungroup()  # Ungroup to return to a standard data frame
  
  
  PSE_not_dups <- PSE_mz1_initial %>% group_by(Metabolite_name) %>% filter(n() == 1)
  
  df <- rbind(PSE_remaining_dups_fixed, PSE_not_dups)
  
  ######pivot longer#####
  df <- df %>% pivot_longer(cols = starts_with('JO_S'), names_to = 'Sample', values_to = 'Area') 
  
  df[,'Sample_Type'] = NA
  df[,'Coating'] = NA
  df[,'Location'] = NA
  
  df <- df %>%   mutate(Coating = ifelse(grepl('_U_', Sample), 'Uncoated', Coating)) %>%
    mutate(Coating = ifelse(grepl('_C_', Sample), 'Coated', Coating)) %>%
    mutate(Location = ifelse(grepl('_O_', Sample), 'Office', Location)) %>%
    mutate(Location = ifelse(grepl('_L_', Sample), 'Lab', Location)) %>%
    mutate(Sample_Type = ifelse(grepl('_Bk_', Sample), 'Blank', Sample_Type)) %>%
    mutate(Sample_Type = ifelse(grepl('_S_', Sample), 'Sample', Sample_Type)) 
  df[,'log_Area_Normalized'] = log10(df$Area)
  df} #starts from blank correction and goes through duplicate removal too but only for coated or uncoated lenses
phylogenic_tree_qsur <- function(df){
  Annotated_data_reorganizeddf <- as.data.frame(df) %>%
    dplyr::rename(Kingdom = kingdom) %>%
    dplyr::rename(Superclass = superclass) %>%
    dplyr::rename(Class = class) %>%
    dplyr::rename(Subclass = subclass)
  
  is_list <- sapply(Annotated_data_reorganizeddf, is.list) 
  
  Annotated_data_reorganizeddf2 <- Annotated_data_reorganizeddf %>%
    mutate_at(.vars = c("Kingdom","Superclass", "Class", "Subclass"),
              .funs = list(~ifelse(.=="", NA, as.character(.))))
  
  data2 <- Annotated_data_reorganizeddf2[complete.cases(Annotated_data_reorganizeddf2[, c("Kingdom", "Superclass", "Class", "Subclass")]), ]
  
  clean_data <- data2 %>% filter(Subclass != 'NULL')
  
  clean_data$Subclass <- gsub(",", "&", clean_data$Subclass)
  clean_data$Subclass <- gsub(" ", "_", clean_data$Subclass)
  
  # Construct a hierarchical Newick tree without numerical leaf nodes
  build_newick <- function(data) {
    kingdoms <- unique(data$Kingdom)
    branches <- sapply(kingdoms, function(kingdom) {
      kingdom_data <- data[data$Kingdom == kingdom, ]
      superclasses <- unique(kingdom_data$Superclass)
      superclass_branches <- sapply(superclasses, function(superclass) {
        superclass_data <- kingdom_data[kingdom_data$Superclass == superclass, ]
        classes <- unique(superclass_data$Class)
        class_branches <- sapply(classes, function(class) {
          class_data <- superclass_data[superclass_data$Class == class, ]
          subclasses <- unique(class_data$Subclass)
          # Remove numerical labels (like "1") when constructing the Newick string
          valid_subclasses <- subclasses[!subclasses %in% c("0.5")]
          paste0("(", paste(valid_subclasses, collapse = ","), "):0.5")
        })
        paste0("(", paste(class_branches, collapse = ","), "):0.5")
      })
      paste0("(", paste(superclass_branches, collapse = ","), "):0.5")
    })
    paste0("(", paste(branches, collapse = ","), ");")
  }
  
  # Generate the Newick string
  newick_string <- build_newick(clean_data)
  
  # Read the tree from the Newick string
  tree <- read.tree(text = newick_string)
  
  # Check for NA branch lengths and handle them
  if (is.null(tree$edge.length)) {
    tree$edge.length <- rep(0.5, nrow(tree$edge))  # Assign default lengths if missing
  } else if (any(is.na(tree$edge.length))) {
    tree$edge.length[is.na(tree$edge.length)] <- 0.5  # Replace NA lengths with default value
  }
  
  # Initialize all edge colors as black
  edge_colors <- rep("black", nrow(tree$edge))
  
  # Assign colors to edges based on subclasses
  for (i in seq_along(tree$tip.label)) {
    tip_label <- tree$tip.label[i]  # Current tip label#
    superclass <- unique(clean_data$Superclass[clean_data$Subclass == tip_label])  # Match superclass
    if (length(superclass) > 0 && any(superclass %in% names(superclass_colors2))) {
      # Find edges leading to this tip and assign the color
      tip_edge <- which(tree$edge[, 2] == i)
      edge_colors[tip_edge] <- superclass_colors2[superclass]
    }
  }
  
  tree$tip.label <- gsub('&',',', tree$tip.label)
  tree$tip.label <- gsub('xxx.*','', tree$tip.label)
  
  # Plot the tree with edge colors
  plot(
    tree,
    type = "fan",                     # Circular fan tree
    show.tip.label = TRUE,            # Show labels at the tips
    tip.color = "black",              # Tip label color
    edge.color = edge_colors,         # Colored branches
    edge.width = 2,                   # Make branches thicker
    cex = 0.7                         # Adjust label size
  )
}
phylogenic_tree_qsur_legend <- function(df){
  Annotated_data_reorganizeddf <- as.data.frame(df) %>%
    dplyr::rename(Kingdom = kingdom) %>%
    dplyr::rename(Superclass = superclass) %>%
    dplyr::rename(Class = class) %>%
    dplyr::rename(Subclass = subclass)
  
  is_list <- sapply(Annotated_data_reorganizeddf, is.list) 
  
  Annotated_data_reorganizeddf2 <- Annotated_data_reorganizeddf %>%
    mutate_at(.vars = c("Kingdom","Superclass", "Class", "Subclass"),
              .funs = list(~ifelse(.=="", NA, as.character(.))))
  
  data2 <- Annotated_data_reorganizeddf2[complete.cases(Annotated_data_reorganizeddf2[, c("Kingdom", "Superclass", "Class", "Subclass")]), ]
  
  clean_data <- data2 %>% filter(Subclass != 'NULL') 
  clean_data$Subclass <- gsub(",", "&", clean_data$Subclass)
  clean_data$Subclass <- gsub(" ", "_", clean_data$Subclass)
  colnames(clean_data)
  #get number of each class
  
  n_distinct(clean_data$Subclass)
  n_distinct(clean_data$Superclass)
  n_distinct(clean_data$Class)
  
  # Construct a hierarchical Newick tree without numerical leaf nodes
  build_newick <- function(data) {
    kingdoms <- unique(data$Kingdom)
    branches <- sapply(kingdoms, function(kingdom) {
      kingdom_data <- data[data$Kingdom == kingdom, ]
      superclasses <- unique(kingdom_data$Superclass)
      superclass_branches <- sapply(superclasses, function(superclass) {
        superclass_data <- kingdom_data[kingdom_data$Superclass == superclass, ]
        classes <- unique(superclass_data$Class)
        class_branches <- sapply(classes, function(class) {
          class_data <- superclass_data[superclass_data$Class == class, ]
          subclasses <- unique(class_data$Subclass)
          # Remove numerical labels (like "1") when constructing the Newick string
          valid_subclasses <- subclasses[!subclasses %in% c("0.5")]
          paste0("(", paste(valid_subclasses, collapse = ","), "):0.5")
        })
        paste0("(", paste(class_branches, collapse = ","), "):0.5")
      })
      paste0("(", paste(superclass_branches, collapse = ","), "):0.5")
    })
    paste0("(", paste(branches, collapse = ","), ");")
  }
  
  # Generate the Newick string
  newick_string <- build_newick(clean_data)
  
  # Read the tree from the Newick string
  tree <- read.tree(text = newick_string)
  
  # Check for NA branch lengths and handle them
  if (is.null(tree$edge.length)) {
    tree$edge.length <- rep(0.5, nrow(tree$edge))  # Assign default lengths if missing
  } else if (any(is.na(tree$edge.length))) {
    tree$edge.length[is.na(tree$edge.length)] <- 0.5  # Replace NA lengths with default value
  }
  
  # Initialize all edge colors as black
  edge_colors <- rep("black", nrow(tree$edge))
  
  # Assign colors to edges based on subclasses
  for (i in seq_along(tree$tip.label)) {
    tip_label <- tree$tip.label[i]  # Current tip label#
    superclass <- unique(clean_data$Superclass[clean_data$Subclass == tip_label])  # Match superclass
    if (length(superclass) > 0 && any(superclass %in% names(superclass_colors2))) {
      # Find edges leading to this tip and assign the color
      tip_edge <- which(tree$edge[, 2] == i)
      edge_colors[tip_edge] <- superclass_colors2[superclass]
    }
  }
  
  tree$tip.label <- gsub('&',',', tree$tip.label)
  tree$tip.label <- gsub('xxx.*','', tree$tip.label)
  
  # Plot the tree with edge colors
  legend("bottom", title = 'Function', legend = names(superclass_colors2), fill = superclass_colors2, cex = 0.7, ncol =2)
  
}

######load and clean non-targeted data files######
#import this way because all columns must be characters, otherwise some information will be erased as NA. this is faster than manually setting all columns as characters#
PSE_pilot_csv <- read_csv("/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/MSDIAL.csv", col_types = cols(.default = "c"))

PSE_clean <- ms_dial_clean(PSE_pilot_csv) %>% #run function on file of interest#
  dplyr::rename(JO_Bk_U_O_1_D = JO1AA) %>% #rename samples so they're more informative, where Bk is blank, S is sample, U is polycarbonate, C is PDMS, O is office, L is lab, _D is to distinguish from solvents
  dplyr::rename(JO_Bk_C_O_1_D = JO2AA) %>%
  dplyr::rename(JO_S_U_O_1_D = JO1A_a) %>%
  dplyr::rename(JO_S_C_O_1_D = JO2A_a) %>%
  dplyr::rename(JO_S_U_O_2_D = JO3A_a) %>%
  dplyr::rename(JO_S_C_O_2_D = JO4AB) %>%
  dplyr::rename(JO_S_U_O_3_D = JO5A) %>%
  dplyr::rename(JO_S_C_O_3_D = JO6A) %>%
  dplyr::rename(JO_Bk_U_L_1_D = JO1BB) %>%
  dplyr::rename(JO_Bk_C_L_1_D = JO2BB) %>%
  dplyr::rename(JO_S_U_L_1_D = JO1B) %>%
  dplyr::rename(JO_S_C_L_1_D = JO2B) %>% 
  dplyr::rename(JO_S_U_L_2_D = JO3B) %>%
  dplyr::rename(JO_S_C_L_2_D = JO4B) %>%
  dplyr::rename(JO_S_U_L_3_D = JO5B) %>%
  dplyr::rename(JO_S_C_L_3_D = JO6B) %>%
  dplyr::rename(JO_MB_IPA = JOIPA) %>%
  dplyr::rename(JO_MB_ACN = JOACN) %>%
  subset(select = -c(JO1A, JO2A, JO3A, JO4A, JO4AA)) #remove samples not used

######import manual integrations#####
#goal: import and clean agilent quant output
SS_manual_csv1 <- read_csv("/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/Surrogates.csv", col_types = cols(.default = "c"), col_names = F) 

SS_manual_csv2 <-SS_manual_csv1 %>% #clean data frame
  subset(select = -c(X1, X2, X5, X6, X7)) %>% #remove column not needed
  row_to_names(2, remove_rows_above = FALSE) #make the second row the row names

manual_list <- list() #make empty list

for(i in seq(3,(ncol(SS_manual_csv2)-2), by = 3)){#loop through all columns by groups of three starting from column 3
  SS_manual_csv2[1,i+1] = SS_manual_csv2[1,i] #make first row of the second column equal to name of the compound 
  SS_manual_csv2[1,i+2] = SS_manual_csv2[1,i] #make first row of the third column equal to name of the compound
  df <-SS_manual_csv2[,c(1:2, i:(i+2))] #make a new df with the first two columns (sample names), and the group of three columns
  df[,'Compound'] <- df[1, 3] #make a new column for the compound, equal to the 1st row of the third column
  df <- df %>%
    mutate(Compound = gsub( " Results", "", Compound)) #remove the Results part of the compound name#
  df <- df[-c(1),] #remove first row now that it's in its own column
  manual_list[[i]] <- df #add all dfs to a list
}

SS_manual_csv <- bind_rows(manual_list) %>% #turn cleaned list into a df
  dplyr::rename(Data_File = `Data File`) %>%
  mutate(Data_File = gsub( ".d", "", Data_File)) %>%
  subset(select = -c(RT, Height, Name)) %>%
  pivot_wider(names_from = Data_File, values_from = Area) %>%#pivot wider so each column is a sample
  dplyr::rename(JO_Bk_U_O_1_D = JO1AA) %>% #rename samples so they're more informative#
  dplyr::rename(JO_Bk_C_O_1_D = JO2AA) %>%
  dplyr::rename(JO_S_U_O_1_D = JO1A_a) %>%
  dplyr::rename(JO_S_C_O_1_D = JO2A_a) %>%
  dplyr::rename(JO_S_U_O_2_D = JO3A_a) %>%
  dplyr::rename(JO_S_C_O_2_D = JO4AB) %>%
  dplyr::rename(JO_S_U_O_3_D = JO5A) %>%
  dplyr::rename(JO_S_C_O_3_D = JO6A) %>%
  dplyr::rename(JO_Bk_U_L_1_D = JO1BB) %>%
  dplyr::rename(JO_Bk_C_L_1_D = JO2BB) %>%
  dplyr::rename(JO_S_U_L_1_D = JO1B) %>%
  dplyr::rename(JO_S_C_L_1_D = JO2B) %>% 
  dplyr::rename(JO_S_U_L_2_D = JO3B) %>%
  dplyr::rename(JO_S_C_L_2_D = JO4B) %>%
  dplyr::rename(JO_S_U_L_3_D = JO5B) %>%
  dplyr::rename(JO_S_C_L_3_D = JO6B) %>%
  dplyr::rename(JO_MB_IPA = JOIPA) %>%
  dplyr::rename(JO_MB_ACN = JOACN) %>%
  subset(select = -c(JO1A, JO2A, JO3A, JO4A, JO4AA)) %>% #remove columns not needed#
  pivot_longer(cols = starts_with('JO'), names_to = 'Sample', values_to = 'Peak_Area') #pivot back now that column names are renamed

######Calculate recoveries######
#goal: calculate surrogate recoveries
SS_manual_csv$Peak_Area <- as.numeric(SS_manual_csv$Peak_Area)
SS_manual_csv_MB <- SS_manual_csv %>%
  filter(Sample == 'JO_MB_ACN') %>% #filter for method blank used as reference
  dplyr::rename(MB_Area = Peak_Area)

SS_manual_csv_rest <- SS_manual_csv %>% #filter for non-method blanks
  filter(Sample != 'JO_MB_IPA' & Sample != 'JO_MB_ACN') %>%
  dplyr::rename(Sample_Area = Peak_Area)

SS_manual_merge <- merge(SS_manual_csv_MB, SS_manual_csv_rest, by = c('Compound')) #match MB peak area to sample peak areas
nrow(SS_manual_csv_rest)-nrow(SS_manual_merge) #sanity check merge - should output 0

SS_manual_merge[,'Recovery'] = SS_manual_merge$Sample_Area*100/SS_manual_merge$MB_Area #calculate recovery for each surrogate of each sample
SS_manual_merge[,'Coating']= NA #empty column for coating
SS_Average_Recovery <- SS_manual_merge %>% #mark in coating based on whether sample name has _U_ or _C_ in it
  mutate(Coating = ifelse(grepl('_U_', Sample.y), 'Uncoated', Coating)) %>%
  mutate(Coating = ifelse(grepl('_C_', Sample.y), 'Coated', Coating)) %>%
  group_by(Compound, Coating) %>%
  get_summary_stats(Recovery, type = 'common') #recovery statistics
SS_Average_Recovery[,'RSD'] = SS_Average_Recovery$sd*100/SS_Average_Recovery$mean #calculate RSD

#avg per compound
SS_manual_csv$Peak_Area <- as.numeric(SS_manual_csv$Peak_Area)
SS_averages <- SS_manual_csv %>% #average peak area for surrogates
  group_by(Compound) %>%
  get_summary_stats(Peak_Area, type = 'common')
SS_averages[,'RSD'] = SS_averages$sd*100/SS_averages$mean #RSD for surrogate peak area

#from earlier peak area stats it looks like MTPP (named 13C-TPP in file) is the most consistent -> normalize to this

PSE_ISTD <- SS_manual_csv %>% 
  pivot_wider(names_from = Compound, values_from = Peak_Area) %>% 
  subset(select = c(Sample, `13C_TPP`)) %>% #extract desired normalization
  dplyr::rename(ISTD = `13C_TPP`)

####merge with targeted istds######
#goal: normalize all peak areas to 13C-TPhP

#pivot PSE_clean longer
PSE_clean_long <- PSE_clean %>% pivot_longer(cols = starts_with('JO'), values_to = 'Uncorrected_Area', names_to = 'Sample') %>% #pivot longer so each row is a feature in a sample  
  mutate(Uncorrected_Area = ifelse(is.na(Uncorrected_Area), 0, Uncorrected_Area)) #change all NAs to 0

PSE_ISTD_Corrections <- merge(PSE_ISTD, PSE_clean_long, by = c('Sample'), all = TRUE)  #merge with istds
#double check that the number of rows in CLEAN_MERGED_ALL_WITH_DUPS_initial_merged is equal CLEAN_MERGED_ALL_WITH_DUPS_initial here#

#sanity check
test1 <- PSE_ISTD_Corrections %>% filter(is.na(Average_Mz)) %>% nrow() #ensure all istds matched a sample - this value should be 0
test1
test2 <- PSE_ISTD_Corrections %>% filter(is.na(ISTD)) %>%nrow() #ensure all samples matched an ISTD. This value should be 0. - one specific sample didn't match things (leave out) - JO_S_C_O_2_D#
test2

#correct for ISTD
PSE_ISTD_Corrections$Uncorrected_Area <-as.numeric(PSE_ISTD_Corrections$Uncorrected_Area) 
PSE_ISTD_Corrections$ISTD <-as.numeric(PSE_ISTD_Corrections$ISTD)
PSE_ISTD_Corrections <- PSE_ISTD_Corrections %>% mutate(ISTD = ifelse(is.na(ISTD), 1, ISTD)) #change things with no ISTD to just 1 so calculations arent affected#

PSE_ISTD_Corrections[,'ISTD_Corrected_Value'] = PSE_ISTD_Corrections$Uncorrected_Area/PSE_ISTD_Corrections$ISTD #normalize

PSE_ISTD_Corrections2 <- PSE_ISTD_Corrections %>% 
  subset(select = -c(Uncorrected_Area, ISTD)) %>% #remove columns no longer needed 
  pivot_wider(values_from = ISTD_Corrected_Value, names_from = Sample)

startingfeatures <- nrow(PSE_ISTD_Corrections2) #starting number of features after normalizing to ISTDß

######blank correct######
#goal: subtract blank values for all features from samples and filter out features with high blanks
Blank_String <- '_Bk_' #define the string unique to the names of blanks
Sample_String <- '_S_' #define the string unique to the names of samples

PSE_ISTD_Corrections[,'Sample_Type'] = NA #empty columns to put in characheristics of the samples
PSE_ISTD_Corrections[,'Coating'] = NA
PSE_ISTD_Corrections[,'Location'] = NA

PSE_Blank_Corrections <- PSE_ISTD_Corrections %>%   mutate(Coating = ifelse(grepl('_U_', Sample), 'Uncoated', Coating)) %>% #add in sample characteristics
  mutate(Coating = ifelse(grepl('_C_', Sample), 'Coated', Coating)) %>%
  mutate(Location = ifelse(grepl('_O_', Sample), 'Office', Location)) %>%
  mutate(Location = ifelse(grepl('_L_', Sample), 'Lab', Location)) %>%
  mutate(Sample_Type = ifelse(grepl(Blank_String, Sample), 'Blank', Sample_Type)) %>%
  mutate(Sample_Type = ifelse(grepl(Sample_String, Sample), 'Sample', Sample_Type)) 

#make separate data frames for blanks and samples to be merged
PSE_blanks <- PSE_Blank_Corrections %>% filter(Sample_Type == 'Blank') %>%  
  subset(select = c(Sample, ISTD, Alignment_ID, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, SMILES, m.z_similarity, MS2_matched, Total_score, ISTD_Corrected_Value, Sample_Type, Coating, Location))
PSE_Samples <- PSE_Blank_Corrections %>% filter(Sample_Type == 'Sample') %>%
  subset(select = c(Sample, ISTD, Alignment_ID, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, SMILES, m.z_similarity, MS2_matched, Total_score, ISTD_Corrected_Value, Sample_Type, Coating, Location))

#merge blanks with samples, matching by feature, location and material. this will be used in later analysis even though it's in the function. 
PSE_Blanks_Matched <- merge(PSE_blanks, PSE_Samples, by = c('Location', 'Coating', 'Alignment_ID', 'Average_Mz', 'Metabolite_name', 'MS.MS_assigned', 'INCHIKEY', 'SMILES', 'm.z_similarity', 'MS2_matched', 'Total_score'), all = T) %>%
  dplyr::rename(Blank_Area = ISTD_Corrected_Value.x) %>% 
  dplyr::rename(Sample_Area = ISTD_Corrected_Value.y) %>%
  subset(select = -c(Sample.x, Location, Coating, ISTD.x, ISTD.y, Sample_Type.y, Sample_Type.x))

nrow(PSE_Blanks_Matched)-nrow(PSE_Samples) #sanity check - output should be 0

PSE_Wide_Blank_Sub <- blank_correction(PSE_Blank_Corrections)

#######replace with LOD#####
#goal: calculate lod and detection frequency, and replace non-detects with LOD

PSE_Wide_Blank_Sub_LODs2 <- PSE_Wide_Blank_Sub %>% 
  lod_and_percent_detection() %>%
  filter(Percent_Detection > 50 | Percent_Detection == 50) #filter for detection frequency

PSE_Wide_Blank_Sub_LODs2$m.z_similarity <- as.numeric(PSE_Wide_Blank_Sub_LODs2$m.z_similarity)
PSE_Wide_Blank_Sub_LODs <-PSE_Wide_Blank_Sub_LODs2 %>% filter(m.z_similarity > 0.95) #apply m/z filter
PSE_Wide_Blank_Sub_LODs$Sample_Avg =apply(PSE_Wide_Blank_Sub_LODs[grep("JO_S", colnames(PSE_Wide_Blank_Sub_LODs), value=TRUE)],1,mean) #get average across samples for later filtering

######remove duplicates - by mz similarity#####
#goal: remove duplicates across all samples, retaining feature with highest m.z similarity
PSE_Wide_Blank_Sub_LODs$m.z_similarity <- as.numeric(PSE_Wide_Blank_Sub_LODs$m.z_similarity)
PSE_dup_int1 <- remove_duplicates(df = PSE_Wide_Blank_Sub_LODs, factor = m.z_similarity)

#####remove duplicates - by sample average#####
#second round of duplicate removal, retaining feature with highest average peak area across samples
PSE_Dups_removed <- remove_duplicates(df = PSE_dup_int1, factor = Sample_Avg)

PSE_Dups_removed$m.z_similarity <- as.numeric(PSE_Dups_removed$m.z_similarity)
PSE_mz1_initial <-PSE_Dups_removed
PSE_remaining_dups <- PSE_mz1_initial %>% group_by(INCHIKEY) %>% filter(n() > 1)
nrow(PSE_remaining_dups) #should output 0

####Look for MSMS matches######
PSE_MSMS <- PSE_Wide_Blank_Sub_LODs %>% filter(MS2_matched == 'TRUE') #filter for MSMS matches
PSE_MSMS$Total_score <- as.numeric(PSE_MSMS$Total_score)
PSE_MSMS_dup_int_1 <- remove_duplicates(df = PSE_MSMS, factor = Total_score) #filter out duplicates, retaining feature with highest total score
PSE_MSMS_dup_removed <- remove_duplicates(df = PSE_MSMS_dup_int_1, factor = Sample_Avg) #for remaining duplicates, retain feeature with highest peak area average across samples

######volcano plots#####
#goal: visualize remaining data (filtering across all samples) as volcano plot
PSE_mz1_longer <- PSE_mz1 %>% pivot_longer(cols = starts_with('JO_S'), names_to = 'Sample', values_to = 'Area') 

PSE_mz1_longer[,'Sample_Type'] = NA
PSE_mz1_longer[,'Coating'] = NA
PSE_mz1_longer[,'Location'] = NA

PSE_mz1_longer <- PSE_mz1_longer %>%   mutate(Coating = ifelse(grepl('_U_', Sample), 'Uncoated', Coating)) %>%
  mutate(Coating = ifelse(grepl('_C_', Sample), 'Coated', Coating)) %>%
  mutate(Location = ifelse(grepl('_O_', Sample), 'Office', Location)) %>%
  mutate(Location = ifelse(grepl('_L_', Sample), 'Lab', Location)) %>%
  mutate(Sample_Type = ifelse(grepl('_Bk_', Sample), 'Blank', Sample_Type)) %>%
  mutate(Sample_Type = ifelse(grepl('_S_', Sample), 'Sample', Sample_Type)) 
PSE_mz1_longer[,'log_Area_Normalized'] = log10(PSE_mz1_longer$Area) #normalize area

##compare abundnaces for coated and uncoated##
PSE_mz1_mwu <- PSE_mz1_longer %>% 
  subset(select = c(Coating, log_Area_Normalized)) %>%
  pivot_wider(names_from = Coating, values_from =log_Area_Normalized)

#by sampler type
PSE_mwu <- PSE_mz1_longer %>%    
  group_by(Metabolite_name, Coating) %>%  #groupings based on PAS sorbent and surrogate class#  
  nest() %>% #turn data into a list
  spread(key = Coating, value = data) %>% #pivot wider so each column is a coating
  mutate(mwu = map2(Coated, Uncoated, ~{t.test(.x$log_Area_Normalized, .y$log_Area_Normalized) %>% #t test test for RSD#   
      tidy()}), Coated = map(Coated, nrow),    Uncoated = map(Uncoated, nrow)  ) %>%  #map data for coated and uncoated to olumns  
  unnest() 
min(PSE_mwu$p.value)
max(PSE_mwu$p.value)

abundance_average <- PSE_mz1_longer %>%
  group_by(Coating, Metabolite_name) %>%
  dplyr::summarise(mean = mean(Area)) %>% #mean area for all features
  pivot_wider(names_from = Coating, values_from = mean) 

PSE_volc_df <- merge(PSE_mwu, abundance_average, by = c('Metabolite_name')) #merge mean areas with volcano data
PSE_volc_df[,'logFC'] = log2(PSE_volc_df$Coated.y/PSE_volc_df$Uncoated.y) #calculate fold change and other values for volcano plot
PSE_volc_df[,'logp']= -log2(PSE_volc_df$p.value) 
PSE_volc_df[,'fold_status'] = NA
log_p_cuotff = -log2(0.05) #to use in plotting (0.05 = p threshold)

PSE_volc_df <- PSE_volc_df %>%
  mutate(fold_status = case_when((logp > log_p_cuotff & logFC > 0) ~ 'Higher in PDMS-Coated', 
                                (logp > log_p_cuotff & logFC < 0) ~ 'Higher in Polycarbonate', TRUE ~ 'Not Significantly Different')) #new column for colour coding in ggplot volcano

sigdiff <- PSE_volc_df %>%  #get statistics for each colour code ctegory
  group_by(fold_status) %>%
  dplyr::summarise(min = min(p.value), max = max(p.value), n=n())

PSE_volc <- PSE_volc_df %>% #plot volcano
  separate_wider_delim(Metabolite_name, "MS2: ", names = c("Discard", "Metabolite_name"), too_few = c('align_end')) %>% #trying to get rid of 'No MS2' 
  separate_wider_delim(Metabolite_name, "low score: ", names = c("Discard2", "Metabolite_name"), too_few = c('align_end')) %>% #trying to get rid of 'No MS2' 
   ggplot(aes(x = logFC, y = logp, colour = fold_status, shape = fold_status)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_classic() +
  geom_hline(yintercept = log_p_cuotff, linetype = 'dashed', colour = 'black', size = 1) +
  labs(x = expression(atop('Log'[2]*' Fold Difference (PDMS-Coated/Polycarbonate)')), y = expression('-Log'[2]*' p value')) +
  scale_colour_manual(values = c('#FBB4AE', '#B3CDE3', 'grey')) +
  theme(legend.position = c(0.75,0.9))+
  scale_y_continuous(limits = c(0, 15))+
  theme(
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 18),
    title = element_text(size = 25), legend.box.background = element_rect(colour = "black")) +# Adjust font size, add box around legend with legend.box.background
  theme(legend.title=element_blank())
PSE_volc
PSE_volc_df[,'FC_status'] = NA

PSE_mwu_location <- PSE_mz1_longer %>%    
  group_by(Metabolite_name, Location) %>%  
  nest() %>%   
  spread(key = Location, value = data) %>%   
  mutate(mwu = map2(Lab, Office, ~{t.test(.x$log_Area_Normalized, .y$log_Area_Normalized) %>% #mwu test for RSD#   
      tidy()}), Lab = map(Lab, nrow),    Office = map(Office, nrow)) %>%   
  unnest() #warnings for ties -> likely because of the LODs, causes normal distribution to be assumed

#by location (same code as above)
abundance_average_location <- PSE_mz1_longer %>%
  group_by(Location, Metabolite_name) %>%
  dplyr::summarise(mean = mean(Area)) %>%
  pivot_wider(names_from = Location, values_from = mean)
PSE_volc_df_location <- merge(PSE_mwu_location, abundance_average_location, by = c('Metabolite_name'))
PSE_volc_df_location[,'logFC'] = log2(PSE_volc_df_location$Lab.y/PSE_volc_df_location$Office.y)
PSE_volc_df_location[,'logp']= -log2(PSE_volc_df_location$p.value) 
PSE_volc_df_location[,'fold_status'] = NA
log_p_cuotff = -log2(0.05)

PSE_volc_df_location <- PSE_volc_df_location %>%
  mutate(fold_status = case_when((logp > log_p_cuotff & logFC > 0) ~ 'Higher in Lab', 
                                 (logp > log_p_cuotff & logFC < 0) ~ 'Higher in Office', TRUE ~ 'Not Significantly Different'))

PSE_volc_location <- PSE_volc_df_location %>% 
  separate_wider_delim(Metabolite_name, "MS2: ", names = c("Discard", "Metabolite_name"), too_few = c('align_end')) %>% #trying to get rid of 'No MS2' 
  separate_wider_delim(Metabolite_name, "low score: ", names = c("Discard2", "Metabolite_name"), too_few = c('align_end')) %>% #trying to get rid of 'No MS2' 
  ggplot(aes(x = logFC, y = logp, colour = fold_status, shape = fold_status)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_classic() +
  scale_x_continuous(limits = c(-16,16))+
  theme(legend.position = c(0.75,0.9))+
  geom_hline(yintercept = log_p_cuotff, linetype = 'dashed', colour = 'black', size = 1) +
   labs(x = expression(atop('Log'[2]*' Fold Difference (Lab/Office)')), y = expression('-Log'[2]*' p value')) +
  scale_colour_manual(values = c('#FBB4AE', '#B3CDE3', 'grey')) +
  theme(
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 18),
    title = element_text(size = 25), legend.box.background = element_rect(colour = "black")) +# Adjust font size, add box around legend with legend.box.background
  theme(legend.title=element_blank())
PSE_volc_location

PSE_volc_grid <- ggarrange(PSE_volc, PSE_volc_location, labels = c('A', 'B'), nrow = 1, font.label=list(color="black",size=35)) #arrnage volcanos into grid

######heatmap#####
#goal: visualize abundance sas a heatmap
PSE_matrix_norm_wide <- PSE_mz1_longer %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_1_D', 'Polycarb Off 1', Sample)) %>% #rename samples so they aren't very long in the figure
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_1_D', 'PDMS Off 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_3_D', 'PDMS Off 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_2_D', 'Polycarb Off 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_2_D', 'PDMS Off 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_3_D', 'Polycarb Off 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_1_D', 'Polycarb Lab 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_1_D', 'PDMS Lab 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_2_D', 'Polycarb Lab 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_2_D', 'PDMS Lab 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_3_D', 'Polycarb Lab 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_3_D', 'PDMS Lab 3', Sample)) %>% 
  subset(select = c('Sample', 'Metabolite_name', 'log_Area_Normalized')) %>% #select only necessary columns for matrix
  pivot_wider(names_from = Sample, values_from = log_Area_Normalized) #pivot wider to turn into matrix

PSE_matrix_norm <- PSE_matrix_norm_wide %>% 
  ungroup()%>%
  subset(select = -c(Metabolite_name)) %>% #remove unneeded column (was needed previously for pivoting)
  data.matrix() #turn into matrix

PSE_pheatmap_norm <- pheatmap(PSE_matrix_norm, show_rownames = FALSE, fontsize_col = 25, fontsize = 18, border_color=NA, scale = 'row') #plot heatmap, normalizing across each row

PSE_heatmap_data <- PSE_mz1_longer %>% #extract the data used for generating heatmap for later use
  subset(select = c('Sample', 'Metabolite_name', 'm.z_similarity','SMILES', 'INCHIKEY', 'log_Area_Normalized', 'LOD', 'MS2_matched')) %>%
  pivot_wider(names_from = Sample, values_from = log_Area_Normalized)

#second heatmap for MSMS matched features
PSE_heatmap_MSMS_data <- PSE_mz1_longer %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_1_D', 'Polycarb Off 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_1_D', 'PDMS Off 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_3_D', 'PDMS Off 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_2_D', 'Polycarb Off 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_O_2_D', 'PDMS Off 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_O_3_D', 'Polycarb Off 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_1_D', 'Polycarb Lab 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_1_D', 'PDMS Lab 1', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_2_D', 'Polycarb Lab 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_2_D', 'PDMS Lab 2', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_U_L_3_D', 'Polycarb Lab 3', Sample)) %>%
  mutate(Sample = ifelse(Sample == 'JO_S_C_L_3_D', 'PDMS Lab 3', Sample)) %>%
  subset(select = c('Sample', 'Metabolite_name', 'log_Area_Normalized', 'MS2_matched')) %>%
  pivot_wider(names_from = Sample, values_from = log_Area_Normalized) %>%
  filter(MS2_matched == 'TRUE') %>%
  filter(Metabolite_name != '[(4R,5S,6R,6aS,7R,10aR,11bR)-5-acetyloxy-6-hydroxy-10a-methoxy-4,7,11b-trimethyl-9-oxo-1,2,3,4a,5,6,6a,7,11,11a-decahydronaphtho[2,1-f][1]benzofuran-4-yl]methyl acetate') %>% #rename features so they are shorter and filter out features with long names
  filter(Metabolite_name !='NCGC00169681-03_C22H26O10_[(1S)-1-(beta-D-Glucopyranosyloxy)-5-hydroxy-1,4a,5,7a-tetrahydrocyclopenta[c]pyran-7-yl]methyl benzoate') %>%
  filter(Metabolite_name != '5-[5-hydroxy-3-(hydroxymethyl)pentyl]-8a-(hydroxymethyl)-5,6-dimethyl-3,4,4a,6,7,8-hexahydronaphthalene-1-carboxylic acid') %>%
  filter(Metabolite_name != '5-(1,2,4a,5-tetramethyl-7-oxo-3,4,8,8a-tetrahydro-2H-naphthalen-1-yl)-3-methylpentanoic acid') %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17alpha-Methyltestosterone', '17-α-Methyltestosterone', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17-beta-Estradiol', '17-β-Estradiol', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'CocamidopropylBetaine', 'Cocamidopropyl Betaine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'Triethanolamine; LC-ESI-QTOF; MS2; CE', 'Triethanolamine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'BIS(2-ETHYLHEXYL)PHTHALATE', 'Bis(2-ethylhexyl) phthalate', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'epsilon-Decalactone', 'ε-Decalactone', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == "Cimifugin 4'-O-beta-D-glucopyranoside", "Cimifugin 4'-O-β-D-glucopyranoside", Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name =="4',5-dihydroxy-7-methoxyisoflavone, 4'-o-beta-d-galactopyranoside", "Prunetin 4'-o-β-d-galactopyranoside", Metabolite_name)) %>%
  ungroup()

PSE_heatmap_MSMS_data_matrix <-   PSE_heatmap_MSMS_data %>%
  subset(select = -c(MS2_matched, Metabolite_name))  %>%
  data.matrix(frame = .) #turn into matrix
rownames(PSE_heatmap_MSMS_data_matrix) <- PSE_heatmap_MSMS_data$Metabolite_name #make matrix rownames the metabolite names

PSE_heatmap_MSMS <- pheatmap(PSE_heatmap_MSMS_data_matrix, show_rownames = TRUE, fontsize_col = 25, fontsize_row = 18, fontsize = 18, border_color=NA, scale = 'row') #fontsize argument is the base font size within the plot

PSE_heatmap_grid <- ggarrange(PSE_pheatmap_norm$gtable, PSE_heatmap_MSMS$gtable, labels = c('A', 'B'), nrow = 1, widths = c(0.7, 1), font.label=list(color="black",size=35)) #arrange heatmaps in grid

######coated vs uncoated fold diff#####
PSE_mz1_longer$Area <- as.numeric(PSE_mz1_longer$Area)

PSE_Fold_Change_Data <- PSE_mz1_longer %>% 
  group_by(SMILES, INCHIKEY, Location, Coating, Metabolite_name, Average_Mz, MS2_matched) %>%
  dplyr::summarise(mean = mean(Area)) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Coating, values_from = mean) #get average per feature per location and coating

PSE_Fold_Change_Data[,'log_fold_diff'] = log10(PSE_Fold_Change_Data$Coated/PSE_Fold_Change_Data$Uncoated) #log10 fold diff for each feature

PSE_Fold_Change_Data$Average_Mz <- as.numeric(PSE_Fold_Change_Data$Average_Mz)

#general fold change plot
PSE_Fold_plot <- PSE_Fold_Change_Data %>%
  ggplot(aes(x = Average_Mz, y = log_fold_diff))+
  geom_point(alpha = 0.7, size = 4, colour = 'grey')+
  theme_classic()+
  geom_hline(yintercept = 0)+
  labs(x = 'm/z', y = 'Log Fold Difference (Coated/Uncoated)')+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ #make plot title in center of plot#
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), #adjust font size#
        axis.text.y = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 18),
        title = element_text(size = 25), strip.text = element_text(size = 18)) #strip.text for plot labels#

PSE_Fold_plot

#lab fold change plot
PSE_Fold_plot_lab <- PSE_Fold_Change_Data %>%
  mutate(MS2_matched = ifelse(MS2_matched == 'TRUE', 'MS² match', 'No MS² match')) %>%
  mutate(MS2_matched = factor(MS2_matched, levels=c('No MS² match', 'MS² match'))) %>% #manually reorder plots#
  filter(Location == 'Lab') %>%
  ggplot(aes(x = Average_Mz, y = log_fold_diff, colour = MS2_matched, shape = MS2_matched))+
  geom_point(alpha = 0.7, size = 4)+
  scale_colour_manual(values = c('black', '#FBB4AE'))+
  theme_classic()+
  geom_hline(yintercept = 0)+
  labs(x = 'm/z', y = 'Log Fold Difference (PDMS-Coated/Polycarbonate)', title = 'Lab')+
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(-5,4))+
  theme(plot.title = element_text(hjust = 0.5))+ #make plot title in center of plot#
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), #adjust font size#
        axis.text.y = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 18),
        title = element_text(size = 25), strip.text = element_text(size = 18)) #strip.text for plot labels#

PSE_Fold_plot_lab

#office fold change plot
PSE_Fold_plot_office <- PSE_Fold_Change_Data %>%
  filter(Location == 'Office') %>%
  ggplot(aes(x = Average_Mz, y = log_fold_diff, colour = MS2_matched, shape = MS2_matched))+
  geom_point(alpha = 0.7, size = 4)+
  scale_colour_manual(values = c('black','#FBB4AE'))+
  theme_classic()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits = c(-5,4))+
  labs(x = 'm/z', y = 'Log Fold Difference (PDMS-Coated/Polycarbonate)', title = 'Office')+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ #make plot title in center of plot#
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), #adjust font size#
        axis.text.y = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 18),
        title = element_text(size = 25), strip.text = element_text(size = 18)) #strip.text for plot labels#

PSE_Fold_plot_office

#get number of higher in coated vs uncoated at each location
status_summary <- PSE_Fold_Change_Data %>%
  mutate(FC_status = ifelse(log_fold_diff >0, 'higher in coated', 'higher in uncoated')) %>%
  group_by(FC_status, Location) %>%
  dplyr::summarise(n = n()) %>% #get n of each
  pivot_wider(names_from = FC_status, values_from = n)
status_summary[,'percent_higher_uncoated'] = status_summary$`higher in uncoated`*100/(status_summary$`higher in coated`+status_summary$`higher in uncoated`) #calculate as percent

#arrange fold change plots in grid
PSE_fold_plot_grid <- ggarrange(PSE_Fold_plot_lab, PSE_Fold_plot_office, labels = c('A', 'B'), nrow = 1, font.label=list(color="black",size=35))

############just for coated or uncoated############
PSE_Blanks_Matched[,'Sample_Type'] = NA
PSE_Blanks_Matched[,'Coating'] = NA
PSE_Blanks_Matched[,'Location'] = NA
PSE_Blanks_Matched <- PSE_Blanks_Matched %>%
  dplyr::rename(Sample = Sample.y) %>%
  mutate(Coating = ifelse(grepl('_U_', Sample), 'Uncoated', Coating)) %>%
  mutate(Coating = ifelse(grepl('_C_', Sample), 'Coated', Coating)) %>%
  mutate(Location = ifelse(grepl('_O_', Sample), 'Office', Location)) %>%
  mutate(Location = ifelse(grepl('_L_', Sample), 'Lab', Location)) %>%
  mutate(Sample_Type = ifelse(grepl('_Bk_', Sample), 'Blank', Sample_Type)) %>%
  mutate(Sample_Type = ifelse(grepl('_S_', Sample), 'Sample', Sample_Type)) 

PSE_Blanks_Matched_C <- PSE_Blanks_Matched %>% 
  filter(Coating == 'Coated') #filter for coated lenses

PSE_Blanks_Matched_U <- PSE_Blanks_Matched %>% 
  filter(Coating == 'Uncoated') #filter for uncoated lenses

#filtering within coated lenses
PSE_C <- redo_analysis(PSE_Blanks_Matched_C) #filtering within coated lenses
PSE_C_wide <- PSE_C %>% subset(select = c(Area, Sample, MS2_matched, Metabolite_name, INCHIKEY, SMILES, Average_Mz, m.z_similarity, Percent_Detection)) %>% pivot_wider(names_from = Sample, values_from = Area) %>%
  group_by(INCHIKEY) %>% #some are ties and dups are still in
  slice(1)
PSE_C_wide[,'Coating'] = 'PDMS-Coated'

#filtering within uncoated lenses
PSE_U <- redo_analysis(PSE_Blanks_Matched_U) #filtering within uncoated lesnes
PSE_U_wide <- PSE_U %>% subset(select = c(Area, Sample, MS2_matched, Metabolite_name, INCHIKEY, SMILES, Average_Mz, m.z_similarity, Percent_Detection)) %>% pivot_wider(names_from = Sample, values_from = Area) %>%
  group_by(INCHIKEY) %>% #some are ties and dups are still in
  slice(1)
PSE_U_wide[,'Coating'] = 'Polycarbonate'

#subsets only with needed information
PSE_U_wide2 <- PSE_U_wide %>% subset(select = c(MS2_matched, Metabolite_name, Coating, INCHIKEY, SMILES, Average_Mz)) 
PSE_C_wide2<- PSE_C_wide %>% subset(select = c(MS2_matched, Metabolite_name, Coating, INCHIKEY, SMILES, Average_Mz))

#bind together filtering from within coated and uncoated
PSE_individual_bind <- rbind(PSE_C_wide2, PSE_U_wide2) 
length(unique(PSE_individual_bind$INCHIKEY)) #get unique compounds

PSE_individual_bind_unique <- PSE_individual_bind %>% #filter so only one of each inchikey exists to speed up models
  group_by(INCHIKEY, SMILES) %>%
  slice(1)

PSE_ven_list <- PSE_individual_bind %>%
  unstack(INCHIKEY ~ Coating) #turn into list where each list is split by coating

redblue = c('#FBB4AE', '#B3CDE3', 'grey') #colours for venn diagram
PSE_ven <- ggvenn(PSE_ven_list, stroke_size = 0.7,fill_alpha = 0.6,count_column = TRUE,
                  fill_color = redblue, set_name_size = 0, text_size = 6, show_percentage = TRUE)+ #venn diagram between coated and uncoated features
  theme(legend.position="none")+
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 18),
    title = element_text(size = 25), strip.text = element_text(size = 18)) # Customize theme
PSE_ven

PSE_individual_bind$Average_Mz <- as.numeric(PSE_individual_bind$Average_Mz)

#density plot overlaying masses captured by lens types
PSE_individual_mass <- PSE_individual_bind %>%
  ggplot(aes(x = Average_Mz, fill = Coating)) +
  theme_classic()+
  geom_density(alpha = 0.25, colour = 'black')+
  labs(x = expression('m/z'), y = 'Density')+  
  scale_fill_brewer(palette = 'Set1')+
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),title = element_text(size = 15), strip.text = element_text(size = 18))+
  theme(legend.position = c(.8,.8)) 
PSE_individual_mass 

#get min and max masses sampled by each sampler type
PSE_individual_massn <- PSE_individual_bind %>%
  group_by(Coating) %>%
  dplyr::summarise(n = n(), min = min(Average_Mz), max(Average_Mz))
PSE_individual_massn$min
PSE_individual_massn$`max(Average_Mz)`

#find exactly which feature is the max and min mass sampled by polycarbonate
poly_mass <- PSE_individual_bind %>% 
  filter(Coating=='Polycarbonate') 

min_poly <- poly_mass %>% filter(Average_Mz == min(poly_mass$Average_Mz))
max_poly <- poly_mass %>% filter(Average_Mz == max(poly_mass$Average_Mz))

#look for phthalates in data
PSE_individual_bind[,'Class'] = NA
phthalates <-PSE_individual_bind %>%
  filter(MS2_matched == 'TRUE') %>%
  mutate(Class = ifelse(grepl('Phthalate', Metabolite_name), 'Phthalate', Class)) %>%
  mutate(Class = ifelse(grepl('phthalate', Metabolite_name), 'Phthalate', Class)) %>%
  mutate(Class = ifelse(grepl('PHTHALATE', Metabolite_name), 'Phthalate', Class)) %>%
  mutate(Class = ifelse(grepl('Phosphate', Metabolite_name), 'Phosphate', Class)) %>%
  mutate(Class = ifelse(grepl('phosphate', Metabolite_name), 'Phosphate', Class)) %>%
  mutate(Class = ifelse(grepl('PHOSPHATE', Metabolite_name), 'Phosphate', Class)) %>%
  filter(!is.na(Class)) %>%
  group_by(Class)

esters <- PSE_individual_bind %>%
  mutate(Class = ifelse(grepl('Ester', Metabolite_name), 'Ester', Class))  %>%
  mutate(Class = ifelse(grepl('ester', Metabolite_name), 'Ester', Class))  %>%
  mutate(Class = ifelse(grepl('ESTER', Metabolite_name), 'Ester', Class))  %>%
  filter(!is.na(Class)) 

#####Query classyfire#######
C_U_bind <- rbind(PSE_C_wide2, PSE_U_wide2)

C_U_bind_unique <- C_U_bind %>%
  group_by(INCHIKEY) %>%
  slice(1) %>%
  ungroup()
MSMS_features1 <- C_U_bind %>% filter(MS2_matched == 'TRUE')
MSMS_featuresn <MSMS_features1 %>% group_by(INCHIKEY) %>%
  dplyr::summarise(n = n())

MSMS_features <- merge(MSMS_featuresn, MSMS_features1, by = c('INCHIKEY'))

PSE_classy_list2<- list() #make empty list

for(i in seq_along(C_U_bind_unique$INCHIKEY)){
  Classification <- get_classification(C_U_bind_unique$INCHIKEY[i])
  # Check if Classification is NULL before proceeding
  if (!is.null(Classification)) {
    Classification_df <- classification(Classification)
    # Add the INCHIKEY as a new column
    Classification_df[,'Compound'] <- C_U_bind_unique$INCHIKEY[i]
    # Store the tibble in the list
    PSE_classy_list2[[i]] <- Classification_df
  } else {
    # print a message or store a default value if NULL
    message("No classification found for: ", C_U_bind_unique$INCHIKEY[i])
  }
}

PSE_classy2 <- bind_rows(PSE_classy_list2) %>%
  subset(select = -c(CHEMONT)) %>%
  pivot_wider(names_from = Level, values_from = Classification) %>%
  subset(select = c(kingdom, superclass, class, subclass, Compound)) %>%
  dplyr::rename(INCHIKEY = Compound)
#write_csv(PSE_classy2,'/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/ClassyfireR.csv')

C_U_bind_classy_merge_missed <- merge(C_U_bind, PSE_classy2, by = c('INCHIKEY'), all = T) %>%
  filter(is.na(superclass)) %>%
  group_by(INCHIKEY) %>%
  slice(1) #obtain these through web based classyfire
#write_csv(C_U_bind_classy_merge_missed,'/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/ClassyfireWeb.csv')

#start fro here to avoid querying classyfire every time
PSE_classy2 <- read_csv('/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/ClassyfireR.csv')

C_U_bind_classy_merge_missed_fixed1 <- read_csv("/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/ClassyfireWeb.csv", col_types = cols(.default = "c")) %>%
  subset(select = -c(`ClassyFy Status`)) %>%
  dplyr::rename(INCHIKEY = InChIKey) %>%
  dplyr::rename(superclass = Superclass) %>%
  dplyr::rename(class = Class) %>%
  dplyr::rename(subclass = Subclass) %>%
  dplyr::rename(kingdom = Kingdom) %>%
  subset(select = c(kingdom, subclass, class, superclass, INCHIKEY))

C_U_classy_all <- rbind(C_U_bind_classy_merge_missed_fixed1, PSE_classy2)
C_U_bind_classy_merge_all <- merge(C_U_bind, C_U_classy_all, by = c('INCHIKEY')) %>%
  filter(!is.na(superclass))

C_U_bind_classy_merge_all_long <- C_U_bind_classy_merge_all %>%
  pivot_longer(cols = c(kingdom, class, subclass, superclass), names_to = 'Level', values_to = 'Classification') 

#####QSUR predictions#####
PSE_Individual_Unique_QSUR <- unique(c(PSE_U_wide$SMILES, PSE_C_wide$SMILES))

PSE_U_wide_subset <- PSE_U_wide %>% 
  subset(select = c(SMILES, Coating, INCHIKEY))
PSE_C_wide_subset <- PSE_C_wide %>%
  subset(select = c(SMILES, Coating, INCHIKEY))
PSE_both <- rbind(PSE_U_wide_subset,PSE_C_wide_subset) 

qsur <- QSUR::qsur_models()
toxp <- QSUR::calculate_toxprints(PSE_Individual_Unique_QSUR)
PSE_Individual_Unique_QSUR_Predictions <- QSUR::predict_all_in_domain(models=qsur, df=toxp)

#write_csv(PSE_Individual_Unique_QSUR_Predictions,'/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/QSUR.csv')

#import from here to avoid running model every time
PSE_Individual_Unique_QSUR_Predictions <- read_csv('/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/QSUR.csv')

PSE_Individual_Unique_QSUR_Predictions_long <-PSE_Individual_Unique_QSUR_Predictions %>%
  pivot_longer(cols = c(2:40), values_to = 'prob', names_to = 'Function') %>%
  dplyr::rename(SMILES=chemical_id) %>%
  filter(!is.na(prob))

PSE_Individual_Unique_QSUR_Predictions_long_merge <- merge(PSE_Individual_Unique_QSUR_Predictions_long,PSE_both, by = c('SMILES'))

PSE_ctxR_both_graph_all2 <- PSE_Individual_Unique_QSUR_Predictions_long_merge %>%
  group_by(Function, Coating) %>%
  mutate(Function = gsub( "_", " ", Function))%>% #remove all underscores
  dplyr::summarise(n=n()) %>%
  mutate(Function = str_to_sentence(Function)) %>%
  mutate(Coating = factor(Coating, levels=c("PDMS-Coated", 'Polycarbonate'))) %>% #reorder so larger bar is on the left
  ggplot(aes(x = reorder(Function,-n), y = n, fill = Coating)) +  ##indicates x and y axes, fill indicates to split bars by subgroup and to create a legend
  geom_bar(stat = "identity", width = 0.7, alpha = 0.4, position = position_dodge2(preserve = 'single')) + ##describes column spacing etc
  labs(x = "Function", y = "Number of Compounds") +    theme_classic()+
  theme(axis.title.x = element_text(size = 20),axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20), axis.title.y=element_text(size=20), legend.text=element_text(size=20),
        legend.title=element_text(size=20))+ 
  scale_fill_brewer(palette = 'Set1')+
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.25)+ # add n on each bar
  theme(plot.margin = unit(c(0, 0.75, 0, 0.75),"inches"))+
  theme(legend.position = c(.85,.85))+ #move legend on graph
  theme(axis.text.x = element_text(angle = 45, hjust=1)) #rotate x axis labels
PSE_ctxR_both_graph_all2

PSE_Individual_Unique_QSUR_Predictions_long_merge_highest <- PSE_Individual_Unique_QSUR_Predictions_long_merge %>%
  group_by(Coating, INCHIKEY, SMILES) %>%
  slice_max(order_by = prob, n = 1)

classy_qsur_merge <- merge(C_U_classy_all, PSE_Individual_Unique_QSUR_Predictions_long_merge_highest, by = c('INCHIKEY')) %>%
  subset(select = -c(superclass)) %>%
  dplyr::rename(superclass = Function) %>%
  filter(!is.na(subclass)) %>%
  mutate(superclass = gsub( "_", " ", superclass))%>% #remove all underscores
  mutate(superclass = str_to_sentence(superclass)) 
  
nrow(classy_qsur_merge)-nrow(PSE_Individual_Unique_QSUR_Predictions_long_merge_highest)

cols <- c('subclass', 'superclass')
classy_qsur_merge$subclass_new = apply(classy_qsur_merge[ ,cols] , 1 , paste , collapse = "xxx" )

classy_qsur_coated <- classy_qsur_merge %>%
  filter(Coating == 'PDMS-Coated') %>%
  subset(select = -c(subclass)) %>%
  dplyr::rename(subclass = subclass_new)

classy_qsur_uncoated <- classy_qsur_merge %>%
  filter(Coating == 'Polycarbonate') %>%
  subset(select = -c(subclass)) %>%
  dplyr::rename(subclass = subclass_new)

classy_qsur_all <- classy_qsur_merge %>%
  subset(select = -c(subclass)) %>%
  dplyr::rename(subclass = subclass_new)

classy_qsur_n <- merge(C_U_classy_all, PSE_Individual_Unique_QSUR_Predictions_long_merge_highest, by = c('INCHIKEY'))  %>%
  filter(!is.na(kingdom) & !is.na(subclass) & !is.na(class) & !is.na(superclass)) #superclass in this case is function

classy_qsur_n_nofilter <- merge(C_U_classy_all, PSE_Individual_Unique_QSUR_Predictions_long_merge_highest, by = c('INCHIKEY'), all = T) %>%
  filter(!is.na(Coating))
classy_qsur_n_coated <- classy_qsur_coated %>%
  filter(!is.na(kingdom) & !is.na(subclass) & !is.na(class) & !is.na(superclass)) #superclass in this case is function
classy_qsur_n_uncoated <- classy_qsur_uncoated %>%
  filter(!is.na(kingdom) & !is.na(subclass) & !is.na(class) & !is.na(superclass)) #superclass in this case is function

length(unique(classy_qsur_n$INCHIKEY))-length(unique(PSE_Individual_Unique_QSUR_Predictions_long_merge_highest$INCHIKEY))

unique(classy_qsur_n$subclass)
unique(classy_qsur_n$superclass)

fills1 <- c(sapply(c('Set1', 'Pastel1'), function(x) brewer.pal(9,x)))
fills2 <- c(sapply(c('Set2', 'Dark2'), function(x) brewer.pal(5,x)))
superclass_colors2 <- setNames(c(fills1, fills2), unique(classy_qsur_merge$superclass))

phylogenic_tree_qsur(classy_qsur_coated)
phylogenic_tree_qsur_legend(classy_qsur_coated)

phylogenic_tree_qsur(classy_qsur_uncoated)
phylogenic_tree_qsur_legend(classy_qsur_uncoated)

phylogenic_tree_qsur(classy_qsur_all)
phylogenic_tree_qsur_legend(classy_qsur_all)


####try normalizing volc#####
PSE_mz1_longer_norm_stats <- PSE_mz1_longer %>%
  group_by(Sample) %>%
  dplyr::summarise(mean = mean(log_Area_Normalized), sd = sd(log_Area_Normalized))
                   
PSE_mz1_longer_norm <- merge(PSE_mz1_longer, PSE_mz1_longer_norm_stats, by = c('Sample'))

PSE_mz1_longer_norm[,'Norm_log_area'] = (PSE_mz1_longer_norm$log_Area_Normalized-PSE_mz1_longer_norm$mean)/PSE_mz1_longer_norm$sd

PSE_t<- PSE_mz1_longer_norm %>%    
  group_by(Metabolite_name, Coating) %>%  #groupings based on PAS sorbent and surrogate class#  
  nest() %>%   
  spread(key = Coating, value = data) %>%   
  mutate(mwu = map2(Coated, Uncoated, ~{t.test(.x$log_Area_Normalized, .y$log_Area_Normalized) %>% 
      tidy()}), Coated = map(Coated, nrow),    Uncoated = map(Uncoated, nrow)  ) %>%   
  unnest() 
min(PSE_t$p.value)
max(PSE_t$p.value)
test <- PSE_t %>% filter(p.value<0.05)

PSE_t2<- PSE_mz1_longer_norm %>%    
  group_by(Metabolite_name, Coating) %>%  #groupings based on PAS sorbent and surrogate class#  
  nest() %>%   
  spread(key = Coating, value = data) %>%   
  mutate(mwu = map2(Coated, Uncoated, ~{t.test(.x$Norm_log_area, .y$Norm_log_area) %>% 
      tidy()}), Coated = map(Coated, nrow),    Uncoated = map(Uncoated, nrow)  ) %>%   
  unnest() 
test <- PSE_t2 %>% filter(p.value<0.05)

#####LogKOA - OPERA###### 

OPERA_all <- read_csv("/Users/anyaguo/Library/CloudStorage/OneDrive-McMasterUniversity/Anya/Eyewear/OPERA.csv", col_types = cols(.default = "c"))

PSE_individual_bind2 <- PSE_individual_bind %>%
  mutate(Metabolite_name = gsub( "_", " ", Metabolite_name))

min(PSE_individual_bind$Average_Mz)
max(PSE_individual_bind$Average_Mz)

OPERA_all$LogKOA_pred <- as.numeric(OPERA_all$LogKOA_pred)

model_OPERA_merge1 <- merge(PSE_individual_bind2, OPERA_all, by = c('Metabolite_name'))
nrow(PSE_individual_bind2)-nrow(model_OPERA_merge1) #still missing some

#re-merges differences in metabolite name can be accounted for
model_OPERA_merge2 <- model_OPERA_merge1 %>% subset(select = c(INCHIKEY, LogKOA_pred)) %>%
  group_by(INCHIKEY) %>%
  slice(1)
model_OPERA_merge <- merge(PSE_individual_bind2, model_OPERA_merge2, by = c('INCHIKEY'))
nrow(model_OPERA_merge)-nrow(PSE_individual_bind2)

model_OPERA_merge_missed <- merge(PSE_individual_bind2, model_OPERA_merge2, by = c('INCHIKEY'), all = T) %>%
  filter(is.na(LogKOA_pred))
length(unique(model_OPERA_merge_missed$INCHIKEY))
length(unique(model_OPERA_merge$INCHIKEY))

model_OPERA_mergen <- model_OPERA_merge %>%
  group_by(Coating) %>%
  dplyr::summarise(n = n())

PSE_individual_bind2n <- PSE_individual_bind2 %>%
  group_by(Coating) %>%
  dplyr::summarise(n = n())

PSE_KOA_dens_OPERA <- model_OPERA_merge %>% 
  mutate(Coating = ifelse(Coating == 'Uncoated', 'Polycarbonate', Coating)) %>%
  mutate(Coating = ifelse(Coating == 'Coated', 'PDMS', Coating)) %>%
  group_by(Coating, Metabolite_name) %>%
  slice(1) %>% #so only unique features
  ggplot(aes(x = LogKOA_pred, fill = Coating)) +
  theme_classic()+
  theme(legend.position = c(0.2, 0.9))+
  geom_density(alpha = 0.25, colour = 'black')+
  labs(x = expression('Log K'[OA]), y = 'Density')+  
  scale_fill_brewer(palette = 'Set1')+
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),title = element_text(size = 15), strip.text = element_text(size = 18))
PSE_KOA_dens_OPERA

model_OPERA_merge_uniques <- model_OPERA_merge %>%
  group_by(INCHIKEY) %>%
  dplyr::summarise(n = n()) %>%
  filter(n == 1)

KOA_range = max(model_OPERA_merge$LogKOA_pred)-min(model_OPERA_merge$LogKOA_pred)

model_OPERA_merge_unique_C <- merge(model_OPERA_merge_uniques, model_OPERA_merge, by = c('INCHIKEY')) %>%
  filter(Coating == 'PDMS-Coated')

model_OPERA_merge_unique_U <- merge(model_OPERA_merge_uniques, model_OPERA_merge, by = c('INCHIKEY')) %>%
  filter(Coating == 'Polycarbonate')

model_OPERA_merge_unique_C$LogKOA_pred <- as.numeric(model_OPERA_merge_unique_C$LogKOA_pred)

PSE_OPERA_C <- model_OPERA_merge_unique_C %>%
  filter(Coating == 'PDMS-Coated') %>%
  ggplot(aes(x = LogKOA_pred)) +
  geom_density(alpha = 0.25, colour = 'black', fill = '#E41A1C')+
  labs(x = expression('Log K'[OA]), y = 'Density')+  
  scale_fill_brewer(palette = 'Set1')+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),title = element_text(size = 15), strip.text = element_text(size = 18))+
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0,0.85))+
  scale_x_continuous(limits = c(0,12))

PSE_OPERA_U <- model_OPERA_merge_unique_U %>%
  filter(Coating == 'Polycarbonate') %>%
  ggplot(aes(x = LogKOA_pred)) +
  geom_density(alpha = 0.25, colour = 'black', fill = '#377EB8')+
  labs(x = expression('Log K'[OA]), y = 'Density')+  
  scale_fill_brewer(palette = 'Set1')+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),title = element_text(size = 15), strip.text = element_text(size = 18))+
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0,0.85))+
  scale_x_continuous(limits = c(0,12))

PSE_OPERA_unique <- merge(model_OPERA_merge_uniques, model_OPERA_merge, by = c('INCHIKEY')) %>% 
  ggplot(aes(x = LogKOA_pred, fill = Coating)) +
  geom_density(alpha = 0.7, colour = 'black')+
  labs(x = expression('Log K'[OA]), y = 'Density')+  
  scale_fill_manual(values = c('#FBB4AE', '#B3CDE3'))+
  theme_classic()+
#  facet_wrap(~Coating, scales='fixed', dir = "v", strip.position = "top",
#             axes = "all", nrow = 2)+
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),title = element_text(size = 15), strip.text = element_text(size = 18))+
  theme(legend.position = c(0.2,0.9))

physchem_update <- ggarrange(PSE_ven, PSE_OPERA_unique, labels = c('A', 'B'), nrow = 1, ncol = 2, font.label=list(color="black",size=26))

model_OPERA_merge_unique_U_above8 <-  model_OPERA_merge_unique_U %>% filter(LogKOA_pred > 8)
model_OPERA_merge_unique_U_below8 <- model_OPERA_merge_unique_U %>% filter(LogKOA_pred < 8)
percent_above8_U <- nrow(model_OPERA_merge_unique_U_above8)*100/nrow(model_OPERA_merge_unique_U)

model_OPERA_merge_unique_C_above8 <-  model_OPERA_merge_unique_C %>% filter(LogKOA_pred > 8)
model_OPERA_merge_unique_C_below8 <- model_OPERA_merge_unique_C %>% filter(LogKOA_pred < 8)
percent_above8_C <- nrow(model_OPERA_merge_unique_C_above8)*100/nrow(model_OPERA_merge_unique_C)

nrow(model_OPERA_merge_unique_U_below8)/nrow(model_OPERA_merge_unique_C_below8)
nrow(model_OPERA_merge_unique_U_below8)-nrow(model_OPERA_merge_unique_C_below8)

nrow(model_OPERA_merge_unique_C_above8)/nrow(model_OPERA_merge_unique_U_above8)
nrow(model_OPERA_merge_unique_C_above8)-nrow(model_OPERA_merge_unique_U_above8)


####uptake capacity#####
#estimate Kpdms-air from tromp et al. equation 2
model_OPERA_merge[,'Kpdms-air'] = model_OPERA_merge$LogKOA_pred*0.778+0.813

sampler_radius = 2.5 #in cm
sampler_height = 0.0102 #in cm

sampler_V = (pi*(sampler_radius^2)*sampler_height)/1000000 #divide by 10^6 to convert to m^3. volume of cylinder
sampler_SA = ((2*pi*sampler_radius*sampler_height) + (2*pi*sampler_radius^2) - (pi*sampler_radius^2))/100 #Surface area of cylinder minus one circle. divide by 100 to convert to dm^2
  
`t25_1-C/Ceq` = 0.75
`t95_1-C/Ceq` = 0.05

avg_walking_speed = 0.5 #m/s. 1.2 m/s is walking speed, so if alternating between walking and sedentary

corrected_R = (0.922*exp(1.489*avg_walking_speed))*sampler_SA

model_OPERA_merge[,'t25_windspeedR'] = ln(`t25_1-C/Ceq`)*(((-sampler_V*(10^model_OPERA_merge$`Kpdms-air`)))/corrected_R)
model_OPERA_merge[,'t95_windspeedR'] = ln(`t95_1-C/Ceq`)*(((-sampler_V*(10^model_OPERA_merge$`Kpdms-air`)))/corrected_R)

model_OPERA_merge_capacities <- model_OPERA_merge %>%
  filter(Coating == 'PDMS-Coated') %>%
  pivot_longer(cols = c('t25_windspeedR', 't95_windspeedR'), names_to = 'type', values_to = 'days')

model_OPERA_merge_capacities[,'log10days'] = log10(model_OPERA_merge_capacities$days)
model_OPERA_merge_capacities_stats <- model_OPERA_merge_capacities %>%
  group_by(type) %>%
  dplyr::summarise(mean = mean(days), min = min(days), max = max(days))

min(model_OPERA_merge_capacities$`Kpdms-air`)
max(model_OPERA_merge_capacities$`Kpdms-air`)
min(model_OPERA_merge$t25_windspeedR)*24*60*60 #in seconds, min t25
max(model_OPERA_merge$t25_windspeedR)/365 #in years, max t25
min(model_OPERA_merge$t95_windspeedR)*24*60 #in minutes
max(model_OPERA_merge$t95_windspeedR)/365 #in years

t25_linear_reg <- model_OPERA_merge_capacities %>%
  filter(type == 't25_windspeedR')  %>%
  lm(data = ., log10days ~ `Kpdms-air`)

t25_coeffs <- summary(t25_linear_reg)$coefficients %>% #get coefficients in data frame
  as.data.frame()

t25_intercept <- t25_coeffs$Estimate[1] #extract intercept
t25_slope <- t25_coeffs$Estimate[2]

t25_intersect = (log10(5)-t25_intercept)/t25_slope #find where t25 line intersects the 10 days


t95_linear_reg <- model_OPERA_merge_capacities %>%
  filter(type == 't95_windspeedR')  %>%
  lm(data = ., log10days ~ `Kpdms-air`)

t95_coeffs <- summary(t95_linear_reg)$coefficients %>% #get coefficients in data frame
  as.data.frame()

t95_intercept <- t95_coeffs$Estimate[1] #extract intercept
t95_slope <- t95_coeffs$Estimate[2]
t95_intersect = (log10(5)-t95_intercept)/t95_slope #find where t95 line intersects the 10 days

model_OPERA_merge_capacities_eq <- model_OPERA_merge_capacities %>%
  filter(`Kpdms-air` < t95_intersect & `Kpdms-air` < t25_intersect) %>%
  group_by(SMILES) %>%
  slice(1) #remove duplicates from t25 and t95 being separate rows
  
model_OPERA_merge_capacities_curvi <- model_OPERA_merge_capacities %>%
  filter(`Kpdms-air` > t95_intersect & `Kpdms-air` < t25_intersect) %>%
  group_by(SMILES) %>%
  slice(1) #remove duplicates from t25 and t95 being separate rows

model_OPERA_merge_capacities_linear <- model_OPERA_merge_capacities %>%
  filter(`Kpdms-air` > t95_intersect & `Kpdms-air` > t25_intersect) %>%
  group_by(SMILES) %>%
  slice(1) #remove duplicates from t25 and t95 being separate rows

(nrow(model_OPERA_merge_capacities)/2)-nrow(model_OPERA_merge_capacities_eq)-nrow(model_OPERA_merge_capacities_curvi)-nrow(model_OPERA_merge_capacities_linear)

percent_eq <- nrow(model_OPERA_merge_capacities_eq)*100/(nrow(model_OPERA_merge_capacities_eq) + nrow(model_OPERA_merge_capacities_curvi) + nrow(model_OPERA_merge_capacities_linear))
percent_curvi <- nrow(model_OPERA_merge_capacities_curvi)*100/(nrow(model_OPERA_merge_capacities_eq) + nrow(model_OPERA_merge_capacities_curvi) + nrow(model_OPERA_merge_capacities_linear))
percent_linear <- nrow(model_OPERA_merge_capacities_linear)*100/(nrow(model_OPERA_merge_capacities_eq) + nrow(model_OPERA_merge_capacities_curvi) + nrow(model_OPERA_merge_capacities_linear))

PSE_pdms_uptake_capacity <- model_OPERA_merge_capacities %>%
  mutate(type = ifelse(type == 't25_windspeedR', 't₂₅', type)) %>%
  mutate(type = ifelse(type == 't95_windspeedR', 't₉₅', type)) %>%
  ggplot(aes(x = `Kpdms-air`, y = log10days, colour = type, shape = type)) +
  geom_point(size = 5, alpha = 0.5)+
  geom_hline(yintercept = log10(5),linetype = 'dashed')+
  geom_vline(xintercept = t25_intersect,linetype = 'dashed')+
  geom_vline(xintercept = t95_intersect,linetype = 'dashed')+
  theme_classic()+
  theme(legend.position = c(0.05, 0.95))+
  scale_colour_manual(values = c('#FBB4AE', '#B3CDE3'))+
  scale_shape_manual(values = c(17, 15))+
  labs(y = expression('Log'[10]*' Days'), x = expression('Log'[10]*' K'[PDMS-Air]))+
  theme(axis.title.x = element_text(size = 23), axis.title.y = element_text(size = 23), 
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18), 
        legend.title = element_text(size = 0), legend.text = element_text(size = 20),
        title = element_text(size = 15), strip.text = element_text(size = 18))

PSE_pdms_uptake_capacity

PSE_pdms_violin <- model_OPERA_merge_capacities %>% 
  mutate(type = ifelse(type == 't25_windspeedR', 't₂₅', type)) %>%
  mutate(type = ifelse(type == 't95_windspeedR', 't₉₅', type)) %>%
  ggplot(aes(y = log10days, x = type, fill = type)) +
  geom_violin(trim=FALSE, alpha = 0.6)+
  scale_fill_manual(values = c('#FBB4AE', '#B3CDE3'))+
  theme_classic()+
  geom_hline(yintercept = log10(5),linetype = 'dashed')+
  theme(legend.position ='none')+
  labs(y = expression('Log'[10]*' Days'), x = '')+
  theme(axis.title.x = element_text(size = 23), axis.title.y = element_text(size = 23), 
        axis.text.x = element_text(size = 23), axis.text.y = element_text(size = 18), 
        legend.title = element_text(size = 0), legend.text = element_text(size = 18),
        title = element_text(size = 15), strip.text = element_text(size = 18))

PSE_pdms_violin

PSE_PDMS_capacity_grid <- ggarrange(PSE_pdms_uptake_capacity, PSE_pdms_violin, labels = c('A', 'B'), nrow = 1, font.label=list(color="black",size=35))

#####bar graphs#####

PSE_Lab_C <- PSE_Blanks_Matched %>%
  filter(Location == 'Lab' & Coating == 'Coated') 

PSE_Lab_U <- PSE_Blanks_Matched %>%
  filter(Location == 'Lab' & Coating == 'Uncoated')

PSE_Off_C <- PSE_Blanks_Matched %>%
  filter(Location == 'Office' & Coating == 'Coated')

PSE_Off_U <- PSE_Blanks_Matched %>%
  filter(Location == 'Office' & Coating == 'Uncoated')

PSE_Lab_C_analyzed <- redo_analysis(PSE_Lab_C) %>%
  group_by(Sample, INCHIKEY) %>% slice(1)

PSE_Lab_U_analyzed <- redo_analysis(PSE_Lab_U)%>%
  group_by(Sample, INCHIKEY) %>% slice(1)

PSE_Off_C_analyzed <- redo_analysis(PSE_Off_C) %>%
  group_by(Sample, INCHIKEY) %>% slice(1)

PSE_Off_U_analyzed <- redo_analysis(PSE_Off_U) %>%
  group_by(Sample, INCHIKEY) %>% slice(1)


PSE_loc_coat_individuals <- rbind(PSE_Lab_C_analyzed, PSE_Lab_U_analyzed, PSE_Off_C_analyzed, PSE_Off_U_analyzed)

MSMS_bar_new_data <- PSE_loc_coat_individuals %>%
  filter(MS2_matched == TRUE) %>%
  filter(Metabolite_name != '[(4R,5S,6R,6aS,7R,10aR,11bR)-5-acetyloxy-6-hydroxy-10a-methoxy-4,7,11b-trimethyl-9-oxo-1,2,3,4a,5,6,6a,7,11,11a-decahydronaphtho[2,1-f][1]benzofuran-4-yl]methyl acetate') %>%
  filter(Metabolite_name != '(8R,9S,10R,13S,14S,17S)-17-hydroxy-10,13,17-trimethyl-2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthren-3-one') %>%
  filter(Metabolite_name != '5-(1,2,4a,5-tetramethyl-7-oxo-3,4,8,8a-tetrahydro-2H-naphthalen-1-yl)-3-methylpentanoic acid') %>%
  filter(Metabolite_name != '5-[5-hydroxy-3-(hydroxymethyl)pentyl]-8a-(hydroxymethyl)-5,6-dimethyl-3,4,4a,6,7,8-hexahydronaphthalene-1-carboxylic acid') %>%
  filter(Metabolite_name != "4',5-dihydroxy-7-methoxyisoflavone, 4'-o-beta-d-galactopyranoside") %>%
  filter(Metabolite_name != 'NCGC00169681-03_C22H26O10_[(1S)-1-(beta-D-Glucopyranosyloxy)-5-hydroxy-1,4a,5,7a-tetrahydrocyclopenta[c]pyran-7-yl]methyl benzoate') %>%
  group_by(Coating, INCHIKEY, Location) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17alpha-Methyltestosterone', '17-α-Methyltestosterone', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17-beta-Estradiol', '17-β-Estradiol', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'CocamidopropylBetaine', 'Cocamidopropyl Betaine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'Triethanolamine; LC-ESI-QTOF; MS2; CE', 'Triethanolamine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'BIS(2-ETHYLHEXYL)PHTHALATE', 'Bis(2-ethylhexyl) phthalate', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'epsilon-Decalactone', 'ε-Decalactone', Metabolite_name)) %>%
  group_by(Metabolite_name, Location, Coating, INCHIKEY) %>% #grouping by metabolite name fine bc filtered out dup metabolite names
  dplyr::summarise(mean_Area = mean(Area), sd_Area = sd(Area), n = n()) %>%
  filter(mean_Area > 0.05) %>%
  ungroup()
MSMS_bar_new_data[,'RSD_Area'] = MSMS_bar_new_data$sd_Area*100/MSMS_bar_new_data$mean_Area  

MSMS_bar_new <- MSMS_bar_new_data %>%
  filter(RSD_Area <20) %>%
  mutate(Coating = ifelse(Coating == 'Coated', 'PDMS-Coated', Coating)) %>%
  mutate(Coating = ifelse(Coating == 'Uncoated', 'Polycarbonate', Coating)) %>%
  ggplot(aes(y = reorder(Metabolite_name, mean_Area), x = mean_Area, fill = Coating))+
  facet_nested_wrap(~Location, dir = "v", strip.position = "top", 
                    axes = "all", remove_labels = "x", scales = 'fixed', nrow = 1) + #split plots#
  theme(legend.position = "right")+ #position legend at right of plot# 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) + ##describes column spacing etc
  labs(x = "", y = "Compound") +    
  theme_classic()+
  geom_errorbarh(aes(xmin = mean_Area-sd_Area, xmax = mean_Area+sd_Area), stat = "identity", position = position_dodge(preserve = "single"), height = 0.7)+
  scale_fill_manual(values = c('#FBB4AE', '#B3CDE3'))+
  theme(plot.margin = unit(c(0, 0.75, 0, 0.75),"inches"))+
  theme(legend.position = c(0.9, 0.85))+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 15),
        title = element_text(size = 25), strip.text = element_text(size = 20), strip.background = element_rect(colour="black", fill="white")) #adjust text size#

MSMS_bar_new

MSMS_bar_new2_data <- PSE_loc_coat_individuals %>%
  filter(MS2_matched == TRUE) %>%
  filter(Metabolite_name != '[(4R,5S,6R,6aS,7R,10aR,11bR)-5-acetyloxy-6-hydroxy-10a-methoxy-4,7,11b-trimethyl-9-oxo-1,2,3,4a,5,6,6a,7,11,11a-decahydronaphtho[2,1-f][1]benzofuran-4-yl]methyl acetate') %>%
  filter(Metabolite_name != '(8R,9S,10R,13S,14S,17S)-17-hydroxy-10,13,17-trimethyl-2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthren-3-one') %>%
  filter(Metabolite_name != '5-(1,2,4a,5-tetramethyl-7-oxo-3,4,8,8a-tetrahydro-2H-naphthalen-1-yl)-3-methylpentanoic acid') %>%
  filter(Metabolite_name != '5-[5-hydroxy-3-(hydroxymethyl)pentyl]-8a-(hydroxymethyl)-5,6-dimethyl-3,4,4a,6,7,8-hexahydronaphthalene-1-carboxylic acid') %>%
  filter(Metabolite_name != "4',5-dihydroxy-7-methoxyisoflavone, 4'-o-beta-d-galactopyranoside") %>%
  filter(Metabolite_name != 'NCGC00169681-03_C22H26O10_[(1S)-1-(beta-D-Glucopyranosyloxy)-5-hydroxy-1,4a,5,7a-tetrahydrocyclopenta[c]pyran-7-yl]methyl benzoate') %>%
  group_by(Coating, Metabolite_name, Location) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17alpha-Methyltestosterone', '17-α-Methyltestosterone', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == '17-beta-Estradiol', '17-β-Estradiol', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'CocamidopropylBetaine', 'Cocamidopropyl Betaine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'Triethanolamine; LC-ESI-QTOF; MS2; CE', 'Triethanolamine', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'BIS(2-ETHYLHEXYL)PHTHALATE', 'Bis(2-ethylhexyl) phthalate', Metabolite_name)) %>%
  mutate(Metabolite_name = ifelse(Metabolite_name == 'epsilon-Decalactone', 'ε-Decalactone', Metabolite_name)) %>%
  dplyr::summarise(mean_Area = mean(Area), sd_Area = sd(Area)) %>%
  filter(mean_Area < 0.05) %>%
  ungroup()
MSMS_bar_new2_data[,'RSD_Area'] = MSMS_bar_new2_data$sd_Area*100/MSMS_bar_new2_data$mean_Area  

MSMS_bar_new2 <- MSMS_bar_new2_data %>% 
  filter(RSD_Area < 20) %>%
  ggplot(aes(y = reorder(Metabolite_name, mean_Area), x = mean_Area, fill = Coating))+
  facet_nested_wrap(~Location, dir = "v", strip.position = "top", 
                    axes = "all", remove_labels = "x", scales = 'fixed', nrow = 1) + #split plots#
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) + ##describes column spacing etc
  labs(x = "Normalized Peak Area", y = "Compound") +    
  theme_classic()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  geom_errorbarh(aes(xmin = mean_Area-sd_Area, xmax = mean_Area+sd_Area), stat = "identity", position = position_dodge(preserve = "single"), height = 0.7)+
  scale_fill_manual(values = c('#FBB4AE', '#B3CDE3'))+
  theme(plot.margin = unit(c(0, 0.75, 0, 0.75),"inches"))+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 20), legend.text = element_text(size = 20),
        title = element_text(size = 25), strip.text = element_text(size = 20), strip.background = element_rect(colour="black", fill="white")) #adjust text size#

MSMS_bar_new2

MSMS_bar_grid <- ggarrange(MSMS_bar_new, MSMS_bar_new2, align = c("hv"), nrow = 2, ncol = 1, font.label=list(color="black",size=26))

