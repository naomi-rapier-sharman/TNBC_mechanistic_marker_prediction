### FPKM_to_TPM.R
# Also saved as Supplementary_File_S11_CellMiner_RNA_seq_in_vitro_logfcs.R

setwd("D:/Documents/Current_Projects/JMN_TNBC_results_new") 
library(tidyverse)

### FUNCTION FROM: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
fpkmToTpm <- function(fpkm){
  ### ORIGINAL VERSION:
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

quotient = function(num, denom) {
  #log(num)/log(denom)
  num/denom
}


necessito_no_undefined = function(columns) {
  This needs to be replace NA with 0.1 and str_replace_all ^0$ with 0.1
}

narrow_book = function(select_columns, full_columns_frame){
  ### Make sure that all requested column names are unique- none repeated.
  keeping_columns <- unique(select_columns) |> print()
  
  #print(as.character(length(select_columns)))
  print(as.character(length(keeping_columns)))
  
  if (length(select_columns) != length(keeping_columns)){
    print("ATTENTION: Some requested column names are duplicate. Keeping only the first occurence of each duplicate column name.")
  }
  if (typeof(keeping_columns) != "character"){
    print("WARNING: select_columns argument must be vector of type character")
  }
  
  naomicols <- colnames(full_columns_frame)
  dudedex = 1
  #print(keeping_columns)
  #print(typeof(keeping_columns))
  for (i in keeping_columns) {
    print(dudedex)
    print(i)
    dudeloco = which(i == naomicols)
    print(paste0("dudeloco = ", as.character(dudeloco)))
    print(paste0("Column in data frame is #", as.character(dudeloco)))
    if (dudedex == 1) {
      #print(full_columns_frame[,dudeloco])
      nano = full_columns_frame[,dudeloco]
      #print(nano)
      nanobook <- nano
    }else {
      nano = full_columns_frame[,dudeloco]
      #print(nano)
      nanobook = cbind(nanobook, nano)
      #print(nanobook)
    }
    nanobook = as_tibble(nanobook)
    #print(nanobook)
    dudedex = dudedex + 1
  }
  #print(nanobook)
  naomi_data1 <- nanobook
  return(naomi_data1)
}

### FUNCTION: undo Cellminer storage of FPKM as logbase2(FPKM + 1)
raw_FPKM = function(x){
  FPKM = (2^(x)) - 1
}


##### CELLMINER YAAAAAAAA

cellminer = read_tsv("CellMiner_RNA-seq_data_TNBC_only.txt") |>
  mutate(Entrez_gene_id = factor(Entrez_gene_id)) %>%
  #mutate(gene_name = factor(gene_name)) %>%
  #mutate(across(.cols = where(is.double), list(undo_log = ~raw_FPKM(.)))) %>% # applies changes, puts new values in new, renamed columns 
  #mutate(across(.cols = where(is.double), list(undo_log = ~(.)))) %>% # makes new, renamed columns without changing the values
  mutate(across(.cols = where(is.double), ~raw_FPKM(.))) %>% # changes column values without making new, renamed columns
  #mutate(across(.cols = where(is.double), ~replace_na(., 0.001))) %>%
  #mutate(across(.cols = where(is.double), ~as.character(.))) %>%
  #mutate(across(.cols = where(is.character), ~str_replace_all(., "^0$", "0.001"))) %>%
  #mutate(across(.cols = where(is.character), ~as.double(.))) %>%
  mutate(across(.cols = where(is.double), ~fpkmToTpm(.))) %>%
  mutate(across(.cols = where(is.factor), ~as.character(.))) %>%
  mutate(across(.cols = where(is.double), ~na_if(., y = 0))) |>
  print()


gene_name <- cellminer |>
  pivot_longer(
    cols = `MDA-MB-231`:`BT-549`,
    names_to = "sample",
    values_to = "TPM"
  ) |>
  print()

desired_genes = read_tsv("desired_genes.txt") |>
  unique()

final_cellminer = left_join(desired_genes, gene_name) |>
  group_by(gene_name, sample) |>
  summarize(
    TPM
  ) |>
  #pivot_wider(
  #  names_from = gene_name,
  #  values_from = TPM
  #) |>
  pivot_wider(
    names_from = sample,
    values_from = TPM
  ) |>
  #mutate(across(.cols = where(is.double), ~quotient(., GAPDH))) |>
  print()






just_phenotype = read_delim("JMN_TNBC_count_matrix_human.tsv", delim = " ") |>
  print()

sample = pull(just_phenotype, Sample)
phenotype = pull(just_phenotype, sample_type.treatment)
print(sample)
print(phenotype)

phen_key = cbind(sample, phenotype) |>
  as_tibble() |>
  mutate(sample = paste0(sample, ".sra")) |>
  print()

rm(just_phenotype)

intro_my = read_tsv("JMN_TNBC_abundance_matrix_human.tsv") %>%
  print()

print(intro_my)



### BRANCH!!!!!!
healthy_my = full_join(phen_key, intro_my) |>
  filter(phenotype == "Healthy_breast_tissue") |>
  mutate(sample = NULL) |>
  mutate(across(.cols = where(is.double), ~mean(.))) %>%
  mutate(sample_type.treatment = NULL) |>
  #group_by(phenotype) |>
  #summarize() |>
  print()

print(tail(colnames(healthy_my)))

healthy_try <- healthy_my[1,] |>
  pivot_longer(
    cols = A1BG:ZZEF1,
    names_to = "gene_name",
    values_to = "primary_healthy_mean"
  ) |>
  mutate(phenotype = NULL) |>
  print()

print(final_cellminer)

fold_cellminer = left_join(final_cellminer, healthy_try) |>
  #mutate(across(.cols = where(is.character), ~factor(.))) |>
  #mutate(across(.cols = where(is.double), ~replace_na(., 0.1))) |>
  #mutate(across(.cols = where(is.double), ~as.character(.))) |>
  #mutate(across(.cols = where(is.character), ~str_replace_all(., "^0$", "0.1"))) |>
  #mutate(across(.cols = where(is.character), ~as.double(.))) |>
  #mutate(across(.cols = where(is.factor), ~as.character(.))) |>
  mutate(across(.cols = where(is.double), ~quotient(., primary_healthy_mean))) |>
  mutate(across(.cols = where(is.double), ~log2(.))) |>
  mutate(across(.cols = where(is.double), ~round(., digits = 2))) |>
  print()

write_tsv(fold_cellminer, "Cellminer_primary_healthy_log2_foldchange.tsv")









#mean-center with scale function - the scale argument should be set to FALSE
mid_my = full_join(phen_key, intro_my) |>
  filter(phenotype == "TN_Breast_Cancer") |>
  #  mutate(across(.cols = where(is.character), ~factor(.))) |>
  #mutate(across(.cols = where(is.double), ~replace_na(., 0.001))) |>
  #  mutate(across(.cols = where(is.double), ~as.character(.))) |>
  #mutate(across(.cols = where(is.character), ~str_replace_all(., "^0$", "0.001"))) |>
  #  mutate(across(.cols = where(is.character), ~as.double(.))) |>
  #  mutate(across(.cols = where(is.factor), ~as.character(.))) |>
  #mutate(across(.cols = where(is.double), ~quotient(., GAPDH))) |>
  print()

print(mid_my)

late_my <- mid_my |>
  #t(mid_my) |>
  mutate(across(.cols = where(is.double), ~mean(.))) %>%
  as_tibble()
late_my <- late_my[1,] |>
  mutate(phenotype = NULL) %>%
  print()

late_my[1,1] <- "my_data_TNBC_average_TPM_quotient"

print(late_my)

my_name <- colnames(late_my)
my_data <- late_my[1,]
print(my_data)



#these_columns <- as.vector(desired_genes[,1])
these_columns <- pull(desired_genes, gene_name)
these_columns <- append(c("sample"), these_columns)

print(these_columns)
print(typeof(these_columns))
#these_columns <- these_columns[[1]]
final_my = narrow_book(these_columns, late_my) |>
  #mutate(across(.cols = where(is.double), ~(.[]))) %>%
  #mutate(across(.cols = where(is.double), ~mean_center(.))) |>
  #as.data.frame() |>
  #select(where(is.double)) |>
  #summarize(across(everything(), mean)) |>
  #pull(MTHFD1L_GAPDHdiv) |>
  print(n =200)



all_these <- bind_rows(final_my, final_cellminer) |>
  print()

rowls <- colnames(all_these) |>
  print()

side_these <- as.data.frame(all_these) |>
  t() |>
  as.data.frame() |>
  print()

write.table(side_these, file = "JMN_TNBC_and_CellMiner_TPM_quotients.tsv", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")



