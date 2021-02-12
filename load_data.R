### DATASETS TO USE ###
version_to_use <- '20q1'
version_url <- "https://figshare.com/articles/dataset/DepMap_20Q1_Public/11791698/3"
sample_info_url <- "https://ndownloader.figshare.com/files/21781221"
gene_effect_url <- "https://ndownloader.figshare.com/files/21521910"
gene_expression_url <- "https://ndownloader.figshare.com/files/21521940"
gecko_gene_effect_url <- "https://ndownloader.figshare.com/files/14246687"

`Dataset Used` <- c(version_url, sample_info_url, 
                    gene_effect_url, gene_expression_url,
                    gecko_gene_effect_url)
names(`Dataset Used`) <- c(version_to_use, "Sample Info", 
                           "Avana DepMap Gene Effect", "Gene Expression",
                           "GeCKO Gene Effect")
### END DATASETS TO USE ###

### FORMAT SAMPLE INFO ###
if(!dir.exists('data')) { dir.create('data') }

if(!file.exists("data/sample_info.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
{
  sample_info <- fread(sample_info_url)
  fwrite(sample_info, "data/sample_info.csv")
} else
{
  sample_info <- fread("data/sample_info.csv")
}

# Fix annotations.
mf <- sample_info %>%
  filter(!is.na(DepMap_ID) & DepMap_ID != '') %>%
  dplyr::select(DepMap_ID, CCLE_name=CCLE_Name, Type=disease, T1=lineage, T2=lineage_subtype, T3=lineage_sub_subtype, age) %>%
  # Custom relabelling of cell lines
  filter(!grepl('MATCHED_NORMAL_TISSUE', CCLE_name)) %>%
  filter(!grepl('ENGINEERED', CCLE_name)) %>%
  filter(!grepl('CHLA57_BONE', CCLE_name)) %>%
  mutate(Type = case_when(
    T1 == 'bone' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>% 
  mutate(Type = case_when(
    T1 == 'central_nervous_system' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  mutate(Type = ifelse(
    grepl('engineer', T1),
    "Immortalized",
    Type
  )) %>%
  mutate(Type = case_when(
    T1 == 'eye' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  mutate(Type = case_when(
    T1 == 'lung' & T2 != "" ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  mutate(Type = case_when(
    T1 == 'soft_tissue' & T2 != "" ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  mutate(Type = ifelse(
    grepl('engineer', Type),
    "Immortalized",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl('lung', Type),
    "Lung",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl('Nsclc', Type),
    "Non-Small Cell Lung",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl('Sclc', Type),
    "Small Cell Lung",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl('Atrt', Type),
    "Rhabdoid",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl('[Rr]habdoid', Type),
    "Rhabdoid",
    Type
  )) %>%
  mutate(Type = ifelse(
    grepl("[Ee]wing",Type), 
    "Ewing", 
    Type)) %>%
  mutate(Type = ifelse(
    CCLE_name == "TASK1_CENTRAL_NERVOUS_SYSTEM", 
    "Ewing", 
    Type)) %>%
  mutate(Type = ifelse(
    T2=='hepatoblastoma', 
    'Hepatoblastoma', 
    Type)) %>%
  mutate(Type = ifelse(
    T2=='PNET', 
    'CNS PNET', 
    Type)) %>%
  mutate(Type = ifelse(
    T2 == "glioma" & (age <= 21 & !is.na(age)), 
    'Pediatric Glioma', 
    Type)) %>%
  mutate(Type = ifelse(
    T2 == "PNET" & (age <= 21 & !is.na(age)), 
    'Pediatric CNS PNET', 
    Type)) %>%
  mutate(Type = ifelse(
    (T2 == "mixed_germ_cell" | T2 == "teratoma") & (age <= 21 & !is.na(age)), 
    'Pediatric Germ Cell', 
    Type)) %>%
  mutate(Type = ifelse(
    grepl('wilms',T2), 
    'Wilms', 
    Type)) %>%
  mutate(Type = ifelse(
    T2 == "MPNST", 
    'MPNST', 
    Type)) %>%
  mutate(Type = ifelse(
    Type == "Leukemia", 
    ifelse(
      T2 == "ALL",
      ifelse(
        T3 == "b_cell",
        "B-ALL",
        ifelse(
          T3 == "t_cell",
          "T-ALL",
          T2)),
      T2),
    Type)) %>%
  mutate(Type = ifelse(
    Type == "Lymphoma", 
    ifelse(
      T2 == "lymphoma_unspecified",
      ifelse(
        T3 == "b_cell",
        "B Cell lymphoma_unspecified",
        ifelse(
          T3 == "t_cell",
          "T Cell lymphoma_unspecified",
          T2)),
      T2),
    Type)) %>%
  mutate(Type = ifelse(
    Type == "non_hodgkin_lymphoma", 
    paste0(T3, " Lymphoma"),
    Type)) %>%
  mutate(Type = stringr::str_to_title(gsub("_"," ",Type))) %>% 
  mutate(Type = ifelse(
    grepl("B-All",Type), 
    "B-ALL", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("T-All",Type), 
    "T-ALL", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Aml",Type), 
    "AML", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Cml",Type), 
    "CML", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Cll",Type), 
    "CLL", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Atl",Type), 
    "ATL", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Nkc",Type), 
    "NKC Leukemia", 
    Type)) %>%
  mutate(Type = ifelse(
    Type == "Dlbcl Lymphoma", 
    "DLBCL",
    Type)) %>%
  mutate(Type = ifelse(
    Type == "T Cell Alcl Lymphoma", 
    "T Cell ALCL",
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Mpnst",Type), 
    "MPNST", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Cns Pnet",Type), 
    gsub("Cns Pnet","CNS PNET",Type), 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Fibrous Histio",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Undifferentiated Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Pleomorphic Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Uterine Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Thyroid Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Round Cell Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Epithelioid Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    grepl("Dermatofibrosarcoma Protuberans",Type), 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = ifelse(
    Type == "Sarcoma", 
    "Sarcoma", 
    Type)) %>%
  mutate(Type = gsub(" Cancer","",Type)) %>%
  mutate(Type = gsub(" Tumor","",Type)) %>%
  distinct()
### END FORMAT SAMPLE INFO ###


### LOAD DATA FUNCTIONS ###


# GENE EFFECT
load_gene_effect <- function(){
  if(file.exists("data/gene_effect.csv.zip") & !file.exists("data/gene_effect.csv"))
  {
    unzip("data/gene_effect.csv.zip", exdir = "data/")
  }
  
  if(!file.exists("data/gene_effect.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_effect <- fread(gene_effect_url)
    fwrite(gene_effect, "data/gene_effect.csv")
  } else
  {
    gene_effect <- fread("data/gene_effect.csv")
  }
  # Reformat as a matrix.
  gene_effect %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_effect$V1
  gene_effect <- as.matrix(gene_effect[,-1])
  DepMap_to_CCLE <- mf$CCLE_name
  names(DepMap_to_CCLE) <- mf$DepMap_ID
  rownames(gene_effect) <- DepMap_to_CCLE[rownames_temp]
  return(gene_effect)
}

# GECKO GENE EFFECT
load_gecko_gene_effect <- function(){
  if(!file.exists("data/gecko_gene_effect.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gecko_gene_effect <- fread(gecko_gene_effect_url)
    fwrite(gecko_gene_effect, "data/gecko_gene_effect.csv")
  } else
  {
    gecko_gene_effect <- fread("data/gecko_gene_effect.csv")
  }
  # Reformat as a matrix.
  rownames_temp <- gecko_gene_effect$V1
  gecko_gene_effect <- as.matrix(gecko_gene_effect[,-1])
  rownames(gecko_gene_effect) <- rownames_temp
  gecko_gene_effect <- t(gecko_gene_effect)
  return(gecko_gene_effect)
}

# GENE EXPRESSION
load_gene_expression <- function(){
  if(file.exists("data/gene_expression.csv.zip") & !file.exists("data/gene_expression.csv"))
  {
    unzip("data/gene_expression.csv.zip", exdir = "data/")
  }
  
  if(!file.exists("data/gene_expression.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_expression <- fread(gene_expression_url)
    fwrite(gene_expression, "data/gene_expression.csv")
  } else
  {
    gene_expression <- fread("data/gene_expression.csv")
  }
  # Reformat as a matrix.
  gene_expression %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_expression$V1
  gene_expression <- as.matrix(gene_expression[,-1])
  rownames(gene_expression) <- rownames_temp
  return(gene_expression)
}


# MSIGDB GENE SETS
load_msigdb <- function(){
  if(!file.exists("data/msigdb.rds")) # If file does not exist, generate from local files. Otherwise load local rds.
  {
    c2 <- read.table(file = "input_files/c2.all.v7.0.symbols.gmt", sep = ",")
    c2_table <- NULL
    for (i in 1:dim(c2)[1])
    {
      c2_data <- unlist(strsplit(c2[i,], "\t"))
      c2_table <- rbind(c2_table, cbind(c2_data[1], c2_data[3:length(c2_data)]))
    }
    colnames(c2_table) <- c("term", "gene")
    msigdb <- as_tibble(c2_table)

    saveRDS(msigdb, "data/msigdb.rds")
  } else
  {
    msigdb <- readRDS("data/msigdb.rds")
  }
  return(msigdb)
}

### END LOAD DATA FUNCTIONS ###