# Argument 1 -> dataframe for thresholds
# Argument 2 -> Bugs_list_joined.tsv
# Argument 3 -> The info table from Strainphlan
# Argument 4 -> The alignment file
# Argument 5 -> Normalized Phylogenetic distance matrix (nGD) for the respective SGB
# Argument 6 -> GGplot output showing the density for each group (same individual / different individual) with the Youden index and 3rd percentile
# Argument 7 -> table output of the SGB with its taxonomy, Alignment and Threshold values
# Argument 8 -> nGD output with T/F or "too few samples to determine" as output

args = base::commandArgs(trailingOnly=TRUE)

project = args[1]

projects=project

base::library("tidyverse")

# Reading the metadata
md <- rio::import(args[1])

# SGB_ID we work with and its info
suppressWarnings(info_sgb <- readr::read_tsv(args[3], col_names = F, show_col_types = F)) #t__SGB4951
alignlen <- readr::read_tsv(args[4], col_names = F, show_col_types=F)
al1 <- alignlen[base::grep(">", alignlen$X1),]
al2 <- alignlen[!alignlen$X1 %in% al1$X1,] %>% 
  base::as.matrix(.) %>% 
  base::as.vector(.)

SGB_ID <- base::gsub("Clade: t__","", info_sgb$X1[1])

# Reading the taxonomy list
sgb_tax <- rio::import(args[2], header = T) %>% 
  dplyr::select(1) %>%
  dplyr::filter(grepl("t__",.$clade_name)) %>% 
  tidyr::separate(., col = clade_name, into = c("Taxonomy", "SGB"), sep = "\\|t__") %>% 
  dplyr::filter(.$SGB %in% SGB_ID)

# Reading the tsv table of pairwise distances
nGD <- readr::read_tsv(file = args[5], col_names = F, show_col_types = F) %>% 
  dplyr::mutate(X1 = base::gsub("_mp4sam","",X1),
                X2 = base::gsub("_mp4sam","",X2))

## create output table with sharing = true or false based on the created threshold in this script
nGD_sharing <- nGD

nGD <- dplyr::left_join(nGD %>% dplyr::select(sampleID_1 = X1, dplyr::everything()),
                        md %>% dplyr::select(sampleID_1 = Sample_ID,
                                             subjectID_1 = Subject_ID,
                                             Dataset_1 = Dataset,
                                             Sample_type_1 = Sample_type,
                                             Donor_1 = Donor_ID))
nGD <- dplyr::left_join(nGD %>% dplyr::select(sampleID_2 = X2, dplyr::everything()),
                        md %>% dplyr::select(sampleID_2 = Sample_ID,
                                             subjectID_2 = Subject_ID,
                                             Dataset_2 = Dataset,
                                             Sample_type_2 = Sample_type,
                                             Donor_2 = Donor_ID)) %>% 
  dplyr::select(sampleID_1, sampleID_2, subjectID_1, subjectID_2, Sample_type_1, Sample_type_2, Donor_1, Donor_2, Dataset_1, Dataset_2, everything()) 

nGD_prep <- nGD %>% dplyr::mutate(same_individual = base::ifelse(.$subjectID_1 == .$subjectID_2, "same_individual", "different_individual"),
                                  same_donor = base::ifelse(.$Donor_1 == .$Donor_2, "same_donor", "different_donor")) %>% 
  dplyr::mutate(potential_related = base::ifelse(.$Sample_type_1 == "Pre_FMT" & .$Sample_type_2 == "Pre_FMT" | .$Sample_type_1 == "Donor" & .$Sample_type_2 == "Donor" |
                                                   .$Sample_type_1 == "Pre_FMT" & .$Sample_type_2 == "Donor" | .$Sample_type_1 == "Donor" & .$Sample_type_2 == "Pre_FMT" |
                                                   .$Sample_type_1 == "Pre_FMT" & .$Sample_type_2 == "Post_FMT" & .$same_donor == "different_donor" |
                                                   .$Sample_type_1 == "Post_FMT" &  .$Sample_type_2 == "Pre_FMT" & .$same_donor == "different_donor" |
                                                   .$Sample_type_1 == "Post_FMT" & .$Sample_type_2 == "Donor" & .$Donor_1 != .$sampleID_2 |
                                                   .$Sample_type_1 == "Donor" & .$Sample_type_2 == "Post_FMT" & .$sampleID_1 != .$Donor_2 |
                                                   .$Sample_type_1 == "Post_FMT" & .$Sample_type_2 == "Post_FMT" & .$same_donor == "different_donor" | 
                                                   .$Sample_type_1 == "Pre_FMT" & .$Sample_type_2 == "Post_FMT" & .$same_individual == "different_individual" |
                                                   .$Sample_type_1 == "Post_FMT" & .$Sample_type_2 == "Pre_FMT" & .$same_individual == "different_individual",
                                                 "unrelated", "pot.related"))
#basically,there is no relatedness IN MY DATASET if:
# The sample comparisons are donor-donor, or pre_fmt-donor, or pre_fmt-pre_FMT 
# The sample comparisons are pre_FMT-post_FMT with different donors
# if post_FMT-donor combinations where the donor of the post_FMT sample is not the same as the donor sample
# two post_FMT samples with different donors

nGD_training <- nGD_prep %>% 
  dplyr::mutate(potential_related = base::factor(potential_related, levels = c("unrelated", "pot.related")))

#numbers of same / different individuals applied in the training
ind_dat <- base::table(nGD_training$potential_related)
print(ind_dat)

if (ind_dat[1] != 0 & ind_dat[2] !=0){
  res_youden <- cutpointr::cutpointr(data = nGD_training, x = X3, class = potential_related, method = cutpointr::maximize_metric, metric = cutpointr::youden)
  sum_youd_cm <- base::as.data.frame(base::summary(res_youden)$confusion_matrix)
  
  if(length(nGD_training[nGD_training$potential_related=="pot.related",]$potential_related) >= 50){
    quantile_pc <- nGD_training %>% dplyr::filter(potential_related == "unrelated") %>% dplyr::pull(X3) %>% stats::quantile(0.05)
    threshold <- base::min(res_youden$optimal_cutpoint, quantile_pc)
    method = base::ifelse(threshold==quantile_pc, "5p", "Youden")
  } else{
    quantile_pc <- nGD_training %>% dplyr::filter(potential_related == "unrelated") %>% dplyr::pull(X3) %>% stats::quantile(0.03)
    threshold <- quantile_pc
    method = "3p"
  }
  
  p_sgb_thres <- ggplot2::ggplot(data = nGD_training) +
    ggplot2::geom_density(ggplot2::aes(x = X3, fill = potential_related), alpha = 0.5 ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = quantile_pc, color = "quantile_pc"), linetype = "dotted", linewidth = 1) +
    ggplot2::theme_minimal() + 
    ggplot2::xlab("StrainPhlAn nGD") + 
    ggplot2::ylab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(legend.title =  ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::scale_color_manual(name = "statistics", values = c(youden = "blue", quantile_pc = "red"))
  
  grDevices::svg(args[6], height=10, width=10)
  ggplot2::theme_set(ggplot2::theme_gray(base_size = 18))
  base::print(p_sgb_thres)
  grDevices::dev.off()
  
  base::prop.table(base::table(nGD_training %>% dplyr::filter(potential_related == "unrelated") %>% dplyr::pull(X3) >= res_youden$optimal_cutpoint))
  
#Create the table
  if (ind_dat[1]>=50){
    table_output <- base::data.frame(SGB_ID = sgb_tax$SGB,
                                     Taxonomy = sgb_tax$Taxonomy,
                                     N_markers = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of markers selected after filtering:", info_sgb$X1),]),
                                     Alignment_length = base::sum(base::nchar(al2))/base::nrow(al1),
                                     N_samples = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of main samples after filtering:", info_sgb$X1),]),
                                     "Avg._fraction_of_gaps" = base::sum(stringr::str_count(al2, "-")) / base::sum(base::nchar(al2)),
                                     Threshold_value = threshold,
                                     Method = method,
                                     Maximum_of_Youden_index = res_youden$sensitivity+res_youden$specificity-1,
                                     FPR = sum_youd_cm$fp/(sum_youd_cm$fp + sum_youd_cm$tn),
                                     FNR = sum_youd_cm$fn/(sum_youd_cm$fn + sum_youd_cm$tp),
                                     N_related_pairs = base::nrow(nGD_training[nGD_training$potential_related=="pot.related",]),
                                     N_unrelated_pairs = base::nrow(nGD_training[nGD_training$potential_related=="unrelated",]))
    nGD_sharing[SGB_ID] <- base::ifelse(nGD_sharing$X3<=threshold, T,F)
  } else{
    table_output <- base::data.frame(SGB_ID = sgb_tax$SGB,
                                     Taxonomy = sgb_tax$Taxonomy,
                                     N_markers = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of markers selected after filtering:", info_sgb$X1),]),
                                     Alignment_length = base::sum(base::nchar(al2))/base::nrow(al1),
                                     N_samples = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of main samples after filtering:", info_sgb$X1),]),
                                     "Avg._fraction_of_gaps" = base::sum(stringr::str_count(al2, "-")) / base::sum(base::nchar(al2)),
                                     Threshold_value = "",
                                     Method = "",
                                     Maximum_of_Youden_index = "",
                                     FPR = "",
                                     FNR = "",
                                     N_related_pairs = base::nrow(nGD_training[nGD_training$potential_related=="pot.related",]),
                                     N_unrelated_pairs = base::nrow(nGD_training[nGD_training$potential_related=="unrelated",]))
    nGD_sharing[SGB_ID] <- "Too few samples to determine"
  }
} else{
  if (ind_dat[1] !=0) {
    p_sgb_thres <- ggplot2::ggplot(data = nGD_training) +
      ggplot2::geom_density(ggplot2::aes(x = X3, fill = potential_related), alpha = 0.5 ) +
      ggplot2::theme_minimal() + 
      ggplot2::xlab("StrainPhlAn nGD") + 
      ggplot2::ylab("") +
      ggplot2::ggtitle("") +
      ggplot2::theme(legend.title =  ggplot2::element_blank(),
                     legend.position = "bottom")
    
    grDevices::svg(args[6], height=10, width=10)
    ggplot2::theme_set(ggplot2::theme_gray(base_size = 18))
    base::print(p_sgb_thres)
    grDevices::dev.off()
  } else {
    p_sgb_thres <- ggplot2::ggplot(data = nGD_training) +
      ggplot2::geom_density(ggplot2::aes(x = X3, fill = potential_related), alpha = 0.5) +
      ggplot2::scale_fill_manual(values = c("pot.related" = "#00BFC4")) +
      ggplot2::theme_minimal() + 
      ggplot2::xlab("StrainPhlAn nGD") + 
      ggplot2::ylab("") +
      ggplot2::ggtitle("") +
      ggplot2::theme(legend.title =  ggplot2::element_blank(),
                     legend.position = "bottom")
    
    grDevices::svg(args[6], height=10, width=10)
    ggplot2::theme_set(ggplot2::theme_gray(base_size = 18))
    base::print(p_sgb_thres)
    grDevices::dev.off()
  }
  
  table_output <- base::data.frame(SGB_ID = sgb_tax$SGB,
                                   Taxonomy = sgb_tax$Taxonomy,
                                   N_markers = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of markers selected after filtering:", info_sgb$X1),]),
                                   Alignment_length = base::sum(base::nchar(al2))/base::nrow(al1),
                                   N_samples = base::gsub("[^[:digit:]]","", info_sgb[base::grep("Number of main samples after filtering:", info_sgb$X1),]),
                                   "Avg._fraction_of_gaps" = base::sum(stringr::str_count(al2, "-")) / base::sum(base::nchar(al2)),
                                   Threshold_value = "",
                                   Method = "",
                                   Maximum_of_Youden_index = "",
                                   FPR = "",
                                   FNR = "",
                                   N_related_pairs = base::nrow(nGD_training[nGD_training$potential_related=="pot.related",]),
                                   N_unrelated_pairs = base::nrow(nGD_training[nGD_training$potential_related=="unrelated",]))
  nGD_sharing[SGB_ID] <- F
}

rio::export(table_output, args[7])
rio::export(nGD_sharing, args[8])