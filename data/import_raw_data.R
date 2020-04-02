library(dplyr)
library(reshape2)

load("~/ownCloud/Phospho-CD4 backup/output/Tsite.rda")
Tsite$Cluster[is.na(Tsite$Clusters)] <- NA


#path <- paste("~/ownCloud/Phospho-CD4/data/PhosphoProteome-Proteome/", "PhosphoDataWithProtAnalysis_20180919.xlsx", sep = "")
path <- "~/ownCloud/Phospho-CD4 backup/manuscript/manuscript_supTable/Final tables/Supplementary Table 1.xlsx"
df_intensity <- readxl::read_excel(path, sheet = 2, skip = 1, na = "NA")
#format biological replicates names
names(df_intensity) <- gsub("_R1_", "_A_", names(df_intensity))
names(df_intensity) <- gsub("_R3_", "_B_", names(df_intensity))
names(df_intensity) <- gsub("_R4_", "_C_", names(df_intensity))
names(df_intensity) <- gsub("_R5_", "_D_", names(df_intensity))

idx_match <- match(Tsite$psiteID, df_intensity$psiteID)
df_merge <- cbind(Tsite[, -which(names(Tsite)=="GeneID")], 
                  df_intensity[idx_match, 
                               c(which(names(df_intensity) %in% c("GeneID", 
                                                                  "pAnova", 
                                                                  "BestFC", 
                                                                  "BestTimePoint", 
                                                                  "ProportioPassStat", 
                                                                  "WarningProteinFC")), 
                                 grep("^Mean Intensity_MV", names(df_intensity)),
                                 grep("^Log2 Norm", names(df_intensity)))])

idx_intensity <- grep("^Mean Intensity_MV", names(df_merge))
names_intensity <- names(df_merge)[idx_intensity]
s<-strsplit(as.character(names_intensity), split="_")

time <- sapply(s, function(x){x[4]})
time <- factor(time, levels = c("NS", "S15", "S30", "S120", "S300", "S600"))
replicate <- sapply(s, function(x){x[2]})
dataset <- sapply(s, function(x){x[3]})

df_cond <- data.frame(name = names_intensity, time = time, replicate = replicate, dataset = dataset)

all_na_pTyr <- apply(df_merge[, as.character(df_cond$name[df_cond$dataset == "pTyr"])], 1, function(x){sum(!is.na(x))==0})
all_na_TiO2 <- apply(df_merge[, as.character(df_cond$name[df_cond$dataset == "TiO2"])], 1, function(x){sum(!is.na(x))==0})
df_merge[df_merge$Residue == "Y" & !all_na_pTyr, as.character(df_cond$name[df_cond$dataset == "TiO2"])] <- NA
df_merge[df_merge$Residue == "S" & !all_na_TiO2, as.character(df_cond$name[df_cond$dataset == "pTyr"])] <- NA
df_merge[df_merge$Residue == "T" & !all_na_TiO2, as.character(df_cond$name[df_cond$dataset == "pTyr"])] <- NA


##################################################################################################
# Keep only corresponding dataset for "S", "T" and "Y" phospho sites

df_merge$Cluster <- factor(as.numeric(df_merge$Cluster))

##################################################################################################
# Update annotations
df_annot <- pannot::get_annotations_uniprot(id = df_merge$Entry,
                                             columns = c("id",
                                                         "keywords",
                                                         "families",
                                                         "go",
                                                         "go(biological_process)",
                                                         "go(molecular_function)",
                                                         "go(cellular_component)"))

idx_match <- match(df_merge$Entry, df_annot$query_id)
df_merge[["GO"]] <- df_annot[["Gene.ontology..GO."]][idx_match]
df_merge[["GO(biological process)"]] <- df_annot[["Gene.ontology..biological.process."]][idx_match]
df_merge[["GO(molecular function)"]] <- df_annot[["Gene.ontology..molecular.function."]][idx_match]
df_merge[["GO(cellular component)"]] <- df_annot[["Gene.ontology..cellular.component."]][idx_match]
df_merge[["Protein.families"]] <- df_annot[["Protein.families"]][idx_match]
df_merge[["Keywords"]] <- df_annot[["Keywords"]][idx_match]


df_merge[["Kinase_reported_mouse_human"]] <- sapply(df_merge[["Kinase_reported_mouse_human"]], 
                                            function(x){s<-strsplit(x, split=";", fixed = TRUE)[[1]]; paste(s[ sapply(s, nchar) > 0 ], collapse = ";")})

save(df_merge, file = "./data/df_merge.rda")
