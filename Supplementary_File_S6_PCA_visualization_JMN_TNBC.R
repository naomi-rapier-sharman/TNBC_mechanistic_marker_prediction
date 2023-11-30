# 2023.07.26_JMN_TNBC_PCA_for_biomarker_direction.R

install.packages("viridis", "dplyr", "tidyverse")

library(tidyverse)
library(grid)
library(viridis)
setwd("E:/Documents/Current_Projects/JMN_TNBC_results_new") 

### Fetch salmon counts
infile = "JMN_TNBC_new_count_matrix_human_samples_tidy.txt"

plain_salmon = read_tsv(infile) %>%
  print()

print(tail(colnames(plain_salmon)))

#Fetch metadata - study of origin, etc.
narrowing_new = read_tsv("JMN_TNBC_metadata.tsv") %>%
  mutate(Samples = str_replace_all(Samples, ".sra", "")) %>%
  dplyr::rename(Sample = Samples) %>%
  print()

plain_salmon = right_join(narrowing_new, plain_salmon, by = "Sample") %>%
  print()

### Take predicted suptype results and add them back in

subtype_key = read_csv("JMN_TNBC_TNBCtype_results_result.csv") %>%
  print()

subtype_key = subtype_key[,1:2] %>%
  print()

with_subtypes = right_join(subtype_key, plain_salmon, by = "Sample") %>%
  mutate(looking = paste0(subtype, "_", sample_type.treatment)) %>%
  mutate(looking = str_replace_all(looking, "NA_TN_Breast_Cancer", "ERlike")) %>%
  mutate(looking = str_replace_all(looking, "NA_Healthy_breast_tissue", "Healthy")) %>%
  mutate(looking = str_replace_all(looking, "_TN_Breast_Cancer", "")) %>%
  mutate(subtype = looking) %>%
  #dplyr::rename(treatment = sample_type.treatment) %>%
  mutate(looking = NULL) %>%
  mutate(sample_type.treatment = NULL) %>%
  #pull(subtype) %>%
  print()

subtype_key = with_subtypes[,1:2] %>%
  #pull(subtype) %>%
  print()

### Prepare tibble for PCA by removing NAs and saving non-numerical columns separately

salmon = with_subtypes %>%
  drop_na() %>%
  #select(where(~ any(. != 0))) %>%
  print()

all_salmon_samples = pull(salmon, Sample) %>%
  print()

all_salmon_disease = pull(salmon, treatment) %>%
  print()

all_salmon_study = pull(salmon, GEO_study) %>%
  print()

all_salmon_subtype = pull(salmon, subtype) %>%
  print()

print(tail(colnames(salmon)))

all_gene_names = colnames(salmon[6:ncol(salmon)]) %>%
  print()

### Make z-scored tibble with only numbers
salmon_numbers = as_tibble(salmon[6:(ncol(salmon))]) %>%
  #select_if(colSums(.) != 0) %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble() %>%
  print()


all_gene_names = colnames(salmon_numbers[1:ncol(salmon_numbers)]) %>%
  print()

#Retrieve metadata, add to z-score matrix
z_salmon = cbind(all_salmon_samples, all_salmon_disease, all_salmon_study, 
                 all_salmon_subtype, salmon_numbers) %>%
  as_tibble() %>%
  dplyr::rename(sample_type.treatment = all_salmon_disease) %>%
  dplyr::rename(Samples = all_salmon_samples) %>%
  print()

### Run PCA on all samples
salmon_groups_all = prcomp(salmon_numbers, center = FALSE, scale = FALSE) %>%
  print()

###This step re-labels the samples with the sample names prcomp choked on (not numeric)

all_salmon_PCA_nums = salmon_groups_all$x %>%
  as_tibble() %>%
  print()

#Put everything together
all_salmon_metadata = cbind(all_salmon_samples, all_salmon_study, 
                            all_salmon_disease, all_salmon_subtype, all_salmon_PCA_nums) %>%
  as_tibble() %>%
  mutate(all_salmon_disease = str_replace_all(all_salmon_disease, "Healthy_breast_tissue", "Healthy Breast Tissue")) %>%
  mutate(all_salmon_disease = str_replace_all(all_salmon_disease, "TN_Breast_Cancer", "Triple Negative Breast Cancer")) %>%
  dplyr::rename(Samples = all_salmon_samples) %>%
  mutate(cb2 = all_salmon_disease) %>%
  mutate(cb2 = str_replace_all(cb2, "Healthy Breast Tissue", "#1f78b4")) %>%
  mutate(cb2 = str_replace_all(cb2, "Triple Negative Breast Cancer", "#fb9a99")) %>%
  #mutate(pubchem = as.factor(pubchem))%>%
  #pull(cb2) %>%
  print()

print(head(colnames(all_salmon_metadata)))

###Graph 

### All colored by disease status
uso <- ggplot(all_salmon_metadata, aes(x=PC1, y=PC2, color = all_salmon_disease)) +
  geom_point() +
  labs(title = "All Samples PCA by Disease Status") +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  labs(color = "Disease Status") +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())

uso

uso +
  scale_colour_viridis_d()


### By GEO study
geo <- ggplot(all_salmon_metadata, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=all_salmon_study)) +
  labs(title = "All Samples PCA by Study of Origin") +
  theme_bw() +
  labs(color = "GEO Accession") +
  #scale_fill_discrete(labels = c("Healthy Breast", "TNBC")) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())

geo

geo +
  scale_colour_viridis_d()

### By TNBC subtype
sub <- ggplot(all_salmon_metadata, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=all_salmon_subtype)) +
  labs(title = "All Samples PCA by Subtype") +
  theme_bw() +
  labs(color = "subtype") +
  #scale_fill_discrete(labels = c("Healthy Breast", "TNBC")) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())

sub

sub +
  scale_colour_viridis_d()