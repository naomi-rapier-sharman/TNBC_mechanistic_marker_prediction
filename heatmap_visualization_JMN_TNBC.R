# heatmap_visualization_JMN_TNBC.R

#install.packages("reshape2")
setwd("E:/Documents/Current_Projects/JMN_TNBC_results_new") 

### Load libraries
library(tidyverse)
library(heatmap3)
library(ggplot2)
#library(reshape2)


### Load Salmon Matrix Data
genematrix0 = read_tsv("JMN_TNBC_new_count_matrix_human_samples_tidy.txt") %>%
#genematrix0 = read_delim("JMN_TNBC_new_count_matrix_human_samples_tidy.txt", delim = " ") %>%
  dplyr::rename(Disease_State = sample_type.treatment) %>%
  print()

#########


### Save metadata
Sample = pull(genematrix0, Sample) %>%
  print()

Disease_State = pull(genematrix0, Disease_State) %>%
  print()

print(tail(colnames(genematrix0)))

all_gene_names = colnames(genematrix0[2:(ncol(genematrix0) -1)]) %>%
  print()

print(tail(all_gene_names))

explore = pull(genematrix0, CENPF) %>%
  print()

### Make z-scored tibble with only numbers

z_genes = as_tibble(genematrix0[2:(ncol(genematrix0) -1)]) %>%
  #mutate(across(all_of(all_gene_names), as.character)) %>%
  #mutate(across(all_of(all_gene_names), ~ str_replace_all(.x, "^0", "1"))) %>%
  #mutate(across(all_of(all_gene_names), as.numeric)) %>%
  #select(where(~ any(. != 0))) %>%
  log2() %>%
  #select_if(colSums(.) != 0) %>%
  #scale() %>%
  as_tibble() %>%
  #pull
  print()

z_genes = cbind(Sample, Disease_State, z_genes) %>%
  as_tibble() %>%
  print()


####### IMPORTANT ###############
### This is the section you change to flip between z-score by sample and z-score by gene.
#### To z-score by sample, # the select line and un# the two t() lines
#### To z-score by gene, # the two t() lines and un# the select line

z_samples = as_tibble(genematrix0[2:(ncol(genematrix0) -1)]) %>%
  mutate(across(all_of(all_gene_names), as.character)) %>%
  mutate(across(all_of(all_gene_names), ~ str_replace_all(.x, "^0", "1"))) %>%
  mutate(across(all_of(all_gene_names), as.numeric)) %>%
  select(where(~ any(. != 0))) %>%
  #t() %>%
  log10() %>%
  scale() %>%
  #t() %>%
  as_tibble() %>%
  print()

z_samples = cbind(Sample, Disease_State, z_samples) %>%
  as_tibble() %>%
  print()






##########
genematrix <- z_samples
#########

print(tail(colnames(genematrix)))

# Get metadata
Sample = genematrix$Sample
Disease_State = genematrix$Disease_State


# Top 10 edgeR genes
ELMOD3 = genematrix$ELMOD3
HJURP = genematrix$HJURP
KIF23 = genematrix$KIF23
KIF11 = genematrix$KIF11
ATAD2 = genematrix$ATAD2
EZH2 = genematrix$EZH2

# Top 5 downregulated biomarkers
CIDEC = genematrix$CIDEC 
CD300LG = genematrix$CD300LG
C14orf180 = genematrix$C14orf180
TNMD = genematrix$TNMD
CFD = genematrix$CFD


# Top 5 upregulated biomarkers - *s represent biomarkers that are also top 10 edgeR
`* ASPM` = genematrix$ASPM
`* RGS1` = genematrix$RGS1
CENPF = genematrix$CENPF
CENPE = genematrix$CENPE
`* KIF14` = genematrix$KIF14


#rm(genematrix)

### Summarize data, keeping only columns I want graphed
### NOTE: Removed FO.1 because it's not in Kennedy's matrix. Also, moved ASPM to the middle of the two groups.
sub_genematrix = cbind(Sample, ELMOD3, HJURP, KIF23, KIF11, ATAD2, EZH2,
                       CIDEC, CD300LG, C14orf180, TNMD, CFD,  
                       `* ASPM`, `* RGS1`, CENPF, CENPE, `* KIF14`, Disease_State) %>%
  as_tibble() %>%
  mutate(Disease_State = str_replace_all(Disease_State, "Healthy_breast_tissue", "Healthy")) %>%
  mutate(Disease_State = str_replace_all(Disease_State, "TN_Breast_Cancer", "TNBC")) %>%
  mutate(ID = paste0(Disease_State, " ", Sample)) %>%
  arrange(ID) %>%
  #mutate(across(sub_genematrix[2:19]),) %>%
  #pull(ID) %>%
  print()

explore = pull(sub_genematrix, Sample) %>%
  print()

row_ID = pull(sub_genematrix, ID) %>%
  print()

numeric_columns = colnames(sub_genematrix) %>%
  print()

numeric_columns = numeric_columns[2:17]

print(numeric_columns)

rev_numeric = rev(numeric_columns) %>%
  print()

#the_dude = melt(sub_genematrix) %>%
#  print()
pre_dude <- sub_genematrix %>%
  print()

#first_half = pre_dude[1:87, ] %>%
#pull(Disease_State) %>%
#  print()

#second_half = pre_dude[88:196, ] %>%
#pull(Disease_State) %>%
#  print()

### NEED TO LEARN PURRR


### Put data in long form for tile graph
the_dude = gather(pre_dude, all_of(numeric_columns), key="Gene", value = "z_center_Expression") %>%
  mutate(Gene = factor(Gene, levels = rev_numeric)) %>%
  mutate(z_center_Expression = as.numeric(z_center_Expression)) %>%
  #pull(ID) %>%
  print()

### Graph heatmap!
ggplot(the_dude, aes(x=ID, y=Gene)) +
  geom_tile(aes(fill = z_center_Expression)) +
  annotate("segment", x=0, xend = 196, y = 10.5, yend=10.5) +
  annotate("segment", x=87, xend = 87, y = .5, yend=16.5) +
  theme(axis.text.x = element_blank()) +
  scale_x_discrete(breaks = NULL) +
  xlab("Healthy                               TNBC") +
  ylab("Top Biomarkers            Top DEGs") +
  coord_cartesian(clip = "off") +
  theme(plot.tag.position = c(.08, .56),
        plot.tag = element_text(size = rel(.75), hjust = 0)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", guide = guide_legend(title = "Expression\n(z-score)")) 
