
### THE FOLLOWING CODE IS ADAPTED FROM Pathway2Targets2.R BY BRETT E. PICKETT

library(tidyverse)
#library(RCurl)
#library(jsonlite)
library(httr)
library(biomaRt)

#install.packages("graphite") THIS IS FROM BIOCONDUCTR - TRY INSTALLING IT FROM THERE
infile = "TNBC_Mechanistic_Markers_to_Check.txt"

#read in file
#remove duplicate ensembl_gene_id
mechanistic_markers = read_tsv(infile) |>
  print()

#get ensembl ids

gene_dude = pull(mechanistic_markers, mechanistic_markers) |>
  print()

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="hgnc_symbol",
  attributes=c("ensembl_gene_id", "hgnc_symbol"),#"definition_1006"),#, "go"),
  values=gene_dude,
  mart=mart)
genes = as_tibble(genes) |>
  dplyr::rename(gene_id = ensembl_gene_id) |>
  dplyr::rename(gene_name = hgnc_symbol) |>
  print()

ensembl_vector = pull(genes, gene_id) |>
  unique() |>
  print()

trac_vector <- as.vector(c(1,2,3,9,10,11,18,19,20,26,27,28))
merged_drugs <- as.data.frame(NULL)
outfile <- paste0(infile,"-Treatments.tsv")
outfile1 <- paste0(infile,"-RankedTargets.tsv")

#define weights
vlow <- 0.01
low <- 0.5
med <- 1
hi <- 2
p3_w <- 1.5
p2_w <- 1
p1_w <- .5
p4_w <- 2



#build query string:
#**note, does not use "tradeNames" or "approvedIndications" as data under "drug" table
query_string = "

       query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      approvedName
      associatedDiseases {
        count
      }
      tractability {
        label
        modality
        value
      }
      safetyLiabilities {
        event
      }
      subcellularLocations{
        location
      }
    knownDrugs {
      uniqueDrugs
      rows {
        drug {
          id
          name
          isApproved
          maximumClinicalTrialPhase
          hasBeenWithdrawn
        }
      }
    }
   }
  }
  "
# Set base URL of GraphQL API endpoint
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

#submit to openTargets.org
for(k in 1:length(ensembl_vector)){
  #k <- 2
  target <- ensembl_vector[k] |> print()
  #target <- "ENSG00000096150" |> print()
  #target <- "ENSG00000139618"  |> print()
  #target<- "ENSG00000141510" |> print()
  #target <- "ENSG00000197919"
  #target <- "ENSG00000185436"
  # Set variables object of arguments to be passed to endpoint
  variables <- list("ensemblId" = target)
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  # Perform POST request
  r <- POST(url=base_url, body=post_body, encode='json')
  # print to console
  
  #print(httr::content(r))
  #print(httr::content(r)$data)
  target_data <- httr::content(r)$data
  #print(target)
  
  #print(target_data[["target"]][["knownDrugs"]][["rows"]])
  #summary(target_data[["target"]][["knownDrugs"]][["rows"]])
  
  #if(length(target_data[["data"]][["drug"]] > 0)){
  if(length(unlist(target_data[["target"]][["knownDrugs"]])) == 0){
    print(paste0("Pathway Member ",k," of ",length(ensembl_vector)," No drugs found: ",ensembl_vector[k]))
    #print("Option 1")
  }else if(length(unlist(target_data[["target"]][["knownDrugs"]])) == 1 & target_data[["target"]][["knownDrugs"]][[1]] == 0){
    print(paste0("Pathway Member ",k," of ",length(ensembl_vector)," No drugs found: ",ensembl_vector[k]))
    #print("Option 2")
  }else{
    print(paste0("Pathway Member ",k," of ",length(ensembl_vector)," Drug(s) found:  ", ensembl_vector[k]))
    
    for(y in 1:length(target_data[["target"]][["knownDrugs"]][["rows"]])){
      y <- 1
      #print(y)
      #print(target_data[["target"]])
      #print(as.vector(unlist(target_data[["target"]][["knownDrugs"]][["rows"]])))
      #print(as.vector(unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[y]][["drug"]])))
      tractability <- 0
      #summary(target_data[["target"]][["tractability"]])
      #print(target_data[["target"]][["tractability"]])
      
      unique_chembl <- as.vector(NULL)
      inter_col <- as.vector(unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[y]][["drug"]]))
      inter_col <- as.character(inter_col)
      #sum tractability at "approved", "advanced clinical trials", or "phase 1 trials" for SM, Ab, PR, and OC
      for(z in 1:length(trac_vector)){
        #z <- 1
        if(length(target_data[["target"]][["tractability"]])==0){
          sm_approved <- "Unknown"
          sm_advancedClinical <- "Unknown"
          sm_phase1 <- "Unknown"
          ab_approved <- "Unknown"
          ab_advancedClinical <- "Unknown"
          ab_phase1 <- "Unknown"
          pr_approved <- "Unknown"
          pr_advancedClinical <- "Unknown"
          pr_phase1 <- "Unknown"
          oc_approved <- "Unknown"
          oc_advancedClinical <- "Unknown"
          oc_phase1 <- "Unknown"
        }else if(length(target_data[["target"]][["tractability"]])>0){
          sm_approved <- target_data[["target"]][["tractability"]][[1]][["value"]]
          sm_advancedClinical <- target_data[["target"]][["tractability"]][[2]][["value"]]
          sm_phase1 <- target_data[["target"]][["tractability"]][[3]][["value"]]
          ab_approved <- target_data[["target"]][["tractability"]][[8]][["value"]]
          ab_advancedClinical <- target_data[["target"]][["tractability"]][[9]][["value"]]
          ab_phase1 <- target_data[["target"]][["tractability"]][[10]][["value"]]
          pr_approved <- target_data[["target"]][["tractability"]][[18]][["value"]]
          pr_advancedClinical <- target_data[["target"]][["tractability"]][[19]][["value"]]
          pr_phase1 <- target_data[["target"]][["tractability"]][[20]][["value"]]
          oc_approved <- target_data[["target"]][["tractability"]][[26]][["value"]]
          oc_advancedClinical <- target_data[["target"]][["tractability"]][[27]][["value"]]
          oc_phase1 <- target_data[["target"]][["tractability"]][[28]][["value"]]
          if(target_data[["target"]][["tractability"]][[z]][["value"]]=="TRUE"){
            tractability <- tractability+1
          }
        }
      }
      if(length(target_data[["target"]][["subcellularLocations"]])==0){
        subCellLoc <- "No data"
      }else if(length(target_data[["target"]][["subcellularLocations"]])>0){
        subCellLoc <- as.vector(unlist(target_data[["target"]][["subcellularLocations"]][[1]]))
      }
      #make sure to only store unique records 
      if(length(unique_chembl) == 0){
        unique_chembl <- inter_col[1]
        # #if drug has no approved indications, then add an N/A for the record
        # if(length(inter_col == 5)){
        #   inter_col[6] <- as.character("N/A")
        # }
        temp_row <- as.vector(c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],target_data[["target"]][["associatedDiseases"]][["count"]],tractability,sm_approved,sm_advancedClinical,sm_phase1,ab_approved,ab_advancedClinical,ab_phase1,pr_approved,pr_advancedClinical,pr_phase1,oc_approved,oc_advancedClinical,oc_phase1,subCellLoc,length(target_data[["target"]][["safetyLiabilities"]]),target_data[["target"]][["knownDrugs"]][["uniqueDrugs"]],inter_col))
        temp_row <- as.data.frame(t(temp_row))
        #print(temp_row)
        
        #if(is.numeric(temp_row[22])==FALSE){
        #  next()
        #}
        merged_drugs <- as.data.frame(rbind(merged_drugs,temp_row),stringsAsFactors = FALSE, drop = FALSE)
        #write.table(df_init, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)
        
      }else if(grep(inter_col[1],unique_chembl) > 0){
        print(paste0("Duplicate drug found for same target: ",inter_col[1]))
        next()
      }else if(grep(inter_col[1],unique_chembl) == 0){
        unique_chembl <- c(unique_chembl, inter_col[1])
        #if drug has no approved indications, then add an N/A for the record
        if(length(inter_col == 5)){
          inter_col[6] <- as.character("N/A")
        }
        temp_row <- as.vector(c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],target_data[["target"]][["associatedDiseases"]][["count"]],tractability,target_data[["target"]][["tractability"]][[1]][["value"]],target_data[["target"]][["tractability"]][[2]][["value"]],target_data[["target"]][["tractability"]][[3]][["value"]],target_data[["target"]][["tractability"]][[9]][["value"]],target_data[["target"]][["tractability"]][[10]][["value"]],target_data[["target"]][["tractability"]][[11]][["value"]],target_data[["target"]][["tractability"]][[18]][["value"]],target_data[["target"]][["tractability"]][[19]][["value"]],target_data[["target"]][["tractability"]][[20]][["value"]],target_data[["target"]][["tractability"]][[26]][["value"]],target_data[["target"]][["tractability"]][[27]][["value"]],target_data[["target"]][["tractability"]][[28]][["value"]],subCellLoc,length(target_data[["target"]][["safetyLiabilities"]]),target_data[["target"]][["knownDrugs"]][["uniqueDrugs"]],inter_col,sig_dbs[i], sig_paths[i]))
        temp_row <- as.data.frame(t(temp_row))
        merged_drugs <- as.data.frame(rbind(merged_drugs,temp_row),stringsAsFactors = FALSE)
      }
      
      rm(inter_col)
      rm(temp_row)
    }
    
    #$temp_df <- as.data.frame(do.call("cbind",list(target_data[["target"]][["id"]], target_data[["target"]][["approvedSymbol"]], target_data[["target"]][["approvedName"]], unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[1]]),sig_dbs[i], sig_paths[i])))
    
    #drug_merged <- as.data.frame(do.call("cbind", list(target_data[["data"]][["target"]][["gene_info"]][["symbol"]],target_data[["data"]][["target"]][["gene_info"]][["name"]],target_data[["data"]][["target"]][["gene_info"]][["geneid"]],target_data[["data"]][["target"]][["activity"]],target_data[["data"]][["disease"]][["name"]],target_data[["data"]][["drug"]][["molecule_name"]],target_data[["data"]][["drug"]][["molecule_type"]],sig_dbs[i],sig_paths[i])))
    ##colnames(drug_merged) <- c("Symbol","Name","Ensembl_ID","Modulation","Disorder","Molecule_Name","Molecule_Type")
    #write.table(drug_merged, file = outfile, row.names = FALSE, col.names=FALSE, sep = "\t", append = TRUE)
  }
  rm(target_data)
}

colnames(merged_drugs) <- c("Target_ID","Target_Symbol","Target_Name","Associated_Disease_Count","Tractability_Count","sm_approved","sm_advancedClinical","sm_phase1","ab_approved","ab_advancedClinical","ab_phase1","pr_approved","pr_advancedClinical","pr_phase1","oc_approved","oc_advancedClinical","oc_phase1","Subcellular_Location","Safety_Liabilities","Number_Unique_Drugs","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn")
print(merged_drugs)
##colnames(merged_drugs) <- c("Target_ID","Target_Symbol","Target_Name","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Pathway_DB","Pathway_Name")
#merged_drugs <- unique(merged_drugs)
#write.table(merged_drugs, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)
write_excel_csv(merged_drugs, outfile)
print("Data gathering...complete")
print("Computing data...")

###########
#sort the existing table by the following order of criteria
# 1)  Pathway_Count,
# 2)  Tractability_Count,
# 2)  Number_Approved,
# 3)  Safety_Liabilities,
# 4)  Number_Unique_Drugs,
# 5)  Associated_Disease_Count,
# 6)  Number_phase3,
# 7)  Number_phase2,
# 8)  Number_phase1,
# 9)  Number_phase4

#Sorting step 1: iterate over whole table
#find counts for each target across all pathways
merged_drugs2 <- subset(merged_drugs, select = -c(Drug_ID,
                                                  Drug_Name,
                                                  Is_FDA_Approved,
                                                  Highest_Clinical_Trial_Phase,
                                                  Has_Been_Withdrawn
))
merged_drugs2 <- unique(merged_drugs2)
#merged_drugs2$Pathway_Name <- as.character(merged_drugs2$Pathway_Name)
merged_drugs2$Target_Symbol <- as.character(merged_drugs2$Target_Symbol)
target1 <- as.vector(merged_drugs2$Target_Symbol)
target1 <- unique(target1)
num_pathways <- as.vector(as.numeric(NULL))
for(a in 1:length(target1)){
  #a <- 3
  merged_drugs_temp <- as.data.frame(NULL)
  merged_drugs_temp <- subset(merged_drugs2, select = c(Target_Symbol))
  #merged_drugs_temp$Pathway_Name <- as.character(merged_drugs_temp$Pathway_Name)
  #merged_drugs_temp$Target_Symbol <- as.character(merged_drugs_temp$Target_Symbol)
  merged_drugs_temp <- subset(merged_drugs_temp, Target_Symbol == target1[a])
  num_pathways[a] <- nrow(merged_drugs_temp)
  #merged_drugs_temp <- subset(merged_drugs2, select = Pathway_Name & Target_Symbol %in% c(target1[a]))
}
rm(merged_drugs_temp)
rm(merged_drugs2)
num_path_df <- as.data.frame(cbind(num_pathways,target1))
num_path_df$num_pathways <- as.numeric(as.character(num_path_df$num_pathways))
num_path_df$target1 <- as.character(num_path_df$target1)
merged_drugs <- merge(merged_drugs,num_path_df, by.x='Target_Symbol', by.y='target1')

#Step2: iterate through each target
merged_drugs3 <- subset(merged_drugs, select = c(Target_Symbol,
                                                 Drug_Name,
                                                 Is_FDA_Approved,
                                                 Highest_Clinical_Trial_Phase,
                                                 Has_Been_Withdrawn
)
)
num4v <- as.vector(as.numeric(NULL))
num3v <- as.vector(as.numeric(NULL))
num2v <- as.vector(as.numeric(NULL))
num1v <- as.vector(as.numeric(NULL))
numApprovedv <- as.vector(as.numeric(NULL))

merged_drugs3 <- unique(merged_drugs3, stringsAsFactors = FALSE)
merged_drugs3$Target_Symbol<- as.character(merged_drugs3$Target_Symbol)
merged_drugs3$Target_Symbol<- as.character(merged_drugs3$Target_Symbol)
merged_drugs3$Is_FDA_Approved<- as.character(merged_drugs3$Is_FDA_Approved)
merged_drugs3$Highest_Clinical_Trial_Phase<- as.numeric(as.character(merged_drugs3$Highest_Clinical_Trial_Phase))

for(b in 1:length(target1)){
  #b <- 1
  phase_temp <- as.data.frame(subset(merged_drugs3, Target_Symbol == target1[b]), stringsAsFactors = FALSE)
  for(c in 1:nrow(phase_temp)){
    numApprovedv[b] <- as.numeric(sum(phase_temp$Is_FDA_Approved == "TRUE"))
    num3v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 3))
    num2v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 2))
    num1v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 1))
    num4v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 4))
  }
}
rm(merged_drugs3)
rm(phase_temp)
num_df <- as.data.frame(cbind(target1,numApprovedv, num3v, num2v, num1v, num4v))
merged_drugs <- merge(merged_drugs,num_df, by.x='Target_Symbol', by.y='target1')

colnames(merged_drugs) <- c("Target_Symbol","Target_ID","Target_Name","Associated_Disease_Count","Tractability_Count","sm_approved","sm_advancedClinical","sm_phase1","ab_approved","ab_advancedClinical","ab_phase1","pr_approved","pr_advancedClinical","pr_phase1","oc_approved","oc_advancedClinical","oc_phase1","Subcellular_Location","Safety_Liabilities","Number_Unique_Drugs","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Target_in_Pathways","num_Approved_Drugs","num_Phase3","num_Phase2","num_Phase1","num_Phase4")
merged_drugs$Tractability_Count <- as.numeric(as.character(merged_drugs$Tractability_Count))
merged_drugs$Safety_Liabilities <- as.numeric(as.character(merged_drugs$Safety_Liabilities))
merged_drugs$Number_Unique_Drugs <- as.numeric(as.character(merged_drugs$Number_Unique_Drugs))
merged_drugs$num_Approved_Drugs <- as.numeric(as.character(merged_drugs$num_Approved_Drugs))
merged_drugs$num_Phase3 <- as.numeric(as.character(merged_drugs$num_Phase3))
merged_drugs$num_Phase2 <- as.numeric(as.character(merged_drugs$num_Phase2))
merged_drugs$num_Phase1 <- as.numeric(as.character(merged_drugs$num_Phase1))
merged_drugs$num_Phase4 <- as.numeric(as.character(merged_drugs$num_Phase4))
merged_drugs$Associated_Disease_Count <- as.numeric(as.character(merged_drugs$Associated_Disease_Count))

#calculate the weighted score by multiplying metrics in each row by weighting factor (low, med, hi)
weighted_score_col <- as.vector(as.numeric(NULL))
for(j in 1:nrow(merged_drugs)){
  #j <- 1
  weighted_score <- as.numeric(0)
  weighted_score <- sum(weighted_score,(merged_drugs$Target_in_Pathways[j]*med))
  weighted_score <- sum(weighted_score,(merged_drugs$Tractability_Count[j]*med))
  weighted_score <- sum(weighted_score,(merged_drugs$num_Approved_Drugs[j]*med))
  weighted_score <- sum(weighted_score,(merged_drugs$Safety_Liabilities[j]*hi*-1))
  weighted_score <- sum(weighted_score,(merged_drugs$Number_Unique_Drugs[j]*med))
  weighted_score <- sum(weighted_score,(merged_drugs$Associated_Disease_Count[j]*med))
  weighted_score <- sum(weighted_score,(merged_drugs$num_Phase3[j]*p3_w*med))
  weighted_score <- sum(weighted_score,(merged_drugs$num_Phase2[j]*p2_w*med))
  weighted_score <- sum(weighted_score,(merged_drugs$num_Phase1[j]*p1_w*med))
  weighted_score <- sum(weighted_score,(merged_drugs$num_Phase4[j]*p4_w*med))
  weighted_score_col[j] <- weighted_score
}
merged_drugs$Weighted_Score <- weighted_score_col
colnames(merged_drugs) <- c("Target_Symbol","Target_ID","Target_Name","Associated_Disease_Count","Tractability_Count","sm_approved","sm_advancedClinical","sm_phase1","ab_approved","ab_advancedClinical","ab_phase1","pr_approved","pr_advancedClinical","pr_phase1","oc_approved","oc_advancedClinical","oc_phase1","Subcellular_Location","Safety_Liabilities","Number_Unique_Drugs","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Target_in_Pathways","num_Approved_Drugs","num_Phase3","num_Phase2","num_Phase1","num_Phase4","Weighted_Score")

merged_drugs <- merged_drugs[order(-merged_drugs$Weighted_Score,
                                   -merged_drugs$Target_in_Pathways,
                                   -merged_drugs$Tractability_Count,
                                   -merged_drugs$num_Approved_Drugs,
                                   merged_drugs$Safety_Liabilities,
                                   -merged_drugs$Number_Unique_Drugs,
                                   -merged_drugs$Associated_Disease_Count,
                                   -merged_drugs$num_Phase3,
                                   -merged_drugs$num_Phase2,
                                   -merged_drugs$num_Phase1,
                                   -merged_drugs$num_Phase4
),]
write.table(merged_drugs, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)

merged_drugs1 <- as.data.frame(merged_drugs, stringsAsFactors = FALSE)
merged_drugs1$Drug_ID <- NULL
merged_drugs1$Is_FDA_Approved <- NULL
merged_drugs1$Highest_Clinical_Trial_Phase <- NULL
merged_drugs1$Has_Been_Withdrawn <- NULL
merged_drugs1$Pathway_DB <- NULL
merged_drugs1$Pathway_Name <- NULL
#merged_drugs1 <- subset(merged_drugs1, select = -c("Drug_ID",
#                                                      "Drug_Name",
#                                                      "Is_FDA_Approved",
#                                                      "Highest_Clinical_Trial_Phase",
#                                                      "Has_Been_Withdrawn",
#                                                      "Pathway_DB",
#                                                      "Pathway_Name"),
#                                      stringsAsFactors = FALSE
#                               )
merged_drugs1$Drug_Name<- NULL
merged_drugs1$Tractability_Count <- as.numeric(as.character(merged_drugs1$Tractability_Count))
merged_drugs1$Safety_Liabilities <- as.numeric(as.character(merged_drugs1$Safety_Liabilities))
merged_drugs1$Number_Unique_Drugs <- as.numeric(as.character(merged_drugs1$Number_Unique_Drugs))
merged_drugs1$num_Approved_Drugs <- as.numeric(as.character(merged_drugs1$num_Approved_Drugs))
merged_drugs1$num_Phase3 <- as.numeric(as.character(merged_drugs1$num_Phase3))
merged_drugs1$num_Phase2 <- as.numeric(as.character(merged_drugs1$num_Phase2))
merged_drugs1$num_Phase1 <- as.numeric(as.character(merged_drugs1$num_Phase1))
merged_drugs1$num_Phase4 <- as.numeric(as.character(merged_drugs1$num_Phase4))
merged_drugs1$Associated_Disease_Count <- as.numeric(as.character(merged_drugs1$Associated_Disease_Count))
merged_drugs1$Weighted_Score <- as.numeric(as.character(merged_drugs1$Weighted_Score))

merged_drugs1 <- unique(merged_drugs1)
merged_drugs1 <- merged_drugs1[order(-merged_drugs1$Weighted_Score,
                                     -merged_drugs1$Target_in_Pathways,
                                     -merged_drugs1$Tractability_Count,
                                     -merged_drugs1$num_Approved_Drugs,
                                     merged_drugs1$Safety_Liabilities,
                                     -merged_drugs1$Number_Unique_Drugs,
                                     -merged_drugs1$Associated_Disease_Count,
                                     -merged_drugs1$num_Phase3,
                                     -merged_drugs1$num_Phase2,
                                     -merged_drugs1$num_Phase1,
                                     -merged_drugs1$num_Phase4
),]
write.table(merged_drugs1, file = outfile1, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)

write_excel_csv(merged_drugs1, paste0(outfile1, ".csv"))

