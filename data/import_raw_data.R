library(dplyr)
library(data.table)

load("~/ownCloud/Phospho-CD4/output/Tsite.rda")
Tsite$Cluster[is.na(Tsite$Clusters)] <- NA
df_intensity <- readxl::read_excel(paste("~/ownCloud/Phospho-CD4/data/PhosphoProteome-Proteome/", 
                                         "PhosphoDataWithProtAnalysis_20180919.xlsx", 
                                         sep = ""), 
                                   sheet = 1, skip = 1, na = "NA")

idx_match <- match(Tsite$psiteID, df_intensity$psiteID)
df_merge <- cbind(Tsite, df_intensity[idx_match, grep("^MeanLoops_", names(df_intensity))])

idx_intensity <- grep("^MeanLoops_", names(df_merge))
names_intensity <- names(df_merge)[idx_intensity]
s<-strsplit(as.character(names_intensity), split="_")

time <- sapply(s, function(x){x[4]})
time <- factor(time, levels = c("NS.", "S15.", "S30.", "S120.", "S300.", "S600."))
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
df_annot <- queryup::get_annotations_uniprot(id = df_merge$Entry,
                                             columns = c("id",
                                                         "keywords",
                                                         "families",
                                                         "go",
                                                         "go(biological_process)",
                                                         "go(molecular_function)",
                                                         "go(cellular_component)"))

idx_match <- match(df_merge$Entry, df_annot$id)
df_merge[["GO terms"]] <- df_annot[["Gene.ontology..GO."]][idx_match]
df_merge[["GO(biological process)"]] <- df_annot[["Gene.ontology..biological.process."]][idx_match]
df_merge[["GO(molecular function)"]] <- df_annot[["Gene.ontology..molecular.function."]][idx_match]
df_merge[["GO(cellular component)"]] <- df_annot[["Gene.ontology..cellular.component."]][idx_match]


save(df_merge, file = "./data/df_merge.rda")