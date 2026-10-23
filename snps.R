##genes from gwascatalog

library(data.table)
library(dplyr)
library(tidyr)

gwas_file <- "data/gwas_catalog_v1.0.2-associations_e114_r2025-08-24.tsv" 
gwas_data <- fread(gwas_file, sep="\t", header=TRUE)

# define hormone-related genes, loop two times
hormone_genes <- c("ESR1", "ESR2", "CYP19A1", "HSD17B1", "AKR1C1")

hormone_genes <- c("HSD17B1", "HSD17B2", "SRD5A2", "SRD5A1", "AR","HSD3B1","HSD3B2")

gwas_long <- gwas_data %>%
  separate_rows(MAPPED_GENE, sep=",\\s*") 

filtered_data <- gwas_long %>%
  filter(MAPPED_GENE %in% hormone_genes)

library(openxlsx)

write.xlsx(filtered_data, "hormone_gene_SNPs.xlsx", rowNames = FALSE)

write.xlsx(filtered_data, "hormone_gene_SNPs1111.xlsx", rowNames = FALSE)

library(readxl)
library(openxlsx)

hormone_data <- read_excel("hormone_gene_SNPsclean.xlsx")

hormone_data <- read_excel("hormone_gene_SNPs1111clean.xlsx")

classify_context <- function(context) {
  if (is.na(context)) return(NA)
  
  if (grepl("missense_variant|stop_gained|splice_donor_5th_base", context, ignore.case = TRUE)) {
    return(1)
  }

  else if (grepl("synonymous_variant|intron_variant|intergenic_variant", context, ignore.case = TRUE)) {
    return(2)
  }

  else if (grepl("3_prime_UTR_variant|regulatory_region_variant|non_coding_transcript_exon_variant", context, ignore.case = TRUE)) {
    return(3)
  }

  else {
    return(NA)
  }
}

hormone_data <- hormone_data %>%
  mutate(Context_class = sapply(CONTEXT, classify_context))

write.xlsx(hormone_data, "hormone_gene_SNPs_with_context.xlsx", rowNames = FALSE)

write.xlsx(hormone_data, "hormone_gene_SNPs_with_context1111.xlsx", rowNames = FALSE)


###UKB genes


library(data.table)
library(SNPassoc)
library(qqman)     
library(ggplot2)   
library(dplyr)
library(bigsnpr)

# read GWAS dataset（.bed, .bim, .fam）
snp_readBed("GWAS/genotype.bed")
snp_obj <- snp_attach("GWAS/genotype.rds")

map <- snp_obj$map
head(map)  


fam <- snp_obj$fam
head(fam) 

G <- snp_obj$genotypes

cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF
allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID  
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])           
colnames(geno_df) <- map$marker.ID     
rownames(geno_df) <- fam$sample.ID    

geno_df$sample.ID <- rownames(geno_df)
geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/SNPs.xlsx")


###AR


library(data.table)
library(SNPassoc)  
library(qqman)     
library(ggplot2)  
library(dplyr)
library(bigsnpr)


snp_readBed("GWAS/ARgenotype.bed")

snp_obj <- snp_attach("GWAS/ARgenotype.rds")

map <- snp_obj$map
head(map)  

fam <- snp_obj$fam
head(fam)  

G <- snp_obj$genotypes


cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

cat("Performing SNP QC...\n")
snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF

allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID 
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])          
colnames(geno_df) <- map$marker.ID   
rownames(geno_df) <- fam$sample.ID    

geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/ARSNPs.xlsx")


#### ESR2


library(data.table)
library(SNPassoc)  
library(qqman)     
library(ggplot2)   
library(dplyr)
library(bigsnpr)


snp_readBed("GWAS/ESR2.bed")

# 载入数据
snp_obj <- snp_attach("GWAS/ESR2.rds")

map <- snp_obj$map
head(map) 


fam <- snp_obj$fam
head(fam)  

G <- snp_obj$genotypes


cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

# SNP QC
cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF
allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID 
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])           
colnames(geno_df) <- map$marker.ID     
rownames(geno_df) <- fam$sample.ID     

geno_df$sample.ID <- rownames(geno_df)
geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/ESR2SNPs.xlsx")

###AKR1C1

library(data.table)
library(SNPassoc)  
library(qqman)     
library(ggplot2)   
library(dplyr)
library(bigsnpr)

snp_readBed("GWAS/AKR1C1.bed")

snp_obj <- snp_attach("GWAS/AKR1C1.rds")

map <- snp_obj$map
head(map)  

fam <- snp_obj$fam
head(fam)  

G <- snp_obj$genotypes


cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

# SNP QC
cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF

allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID  
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])           
colnames(geno_df) <- map$marker.ID      
rownames(geno_df) <- fam$sample.ID     

geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/AKR1C1SNPs.xlsx")


### CYP19A1

library(data.table)
library(SNPassoc)  
library(qqman)    
library(ggplot2)   
library(dplyr)
library(bigsnpr)

snp_readBed("GWAS/CYP19A1.bed")

snp_obj <- snp_attach("GWAS/CYP19A1.rds")

map <- snp_obj$map
head(map)  

fam <- snp_obj$fam
head(fam) 

G <- snp_obj$genotypes

cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

# SNP QC
cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF
allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID 
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])           
colnames(geno_df) <- map$marker.ID      
rownames(geno_df) <- fam$sample.ID      

geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/CYP19A1SNPs.xlsx")

### HSD17B2

library(data.table)
library(SNPassoc)  
library(qqman)    
library(ggplot2)  
library(dplyr)
library(bigsnpr)

snp_readBed("GWAS/HSD17B2.bed")

snp_obj <- snp_attach("GWAS/HSD17B2.rds")

map <- snp_obj$map
head(map) 

fam <- snp_obj$fam
head(fam) 

G <- snp_obj$genotypes

cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

# SNP QC
cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF

allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID  
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])            
colnames(geno_df) <- map$marker.ID      
rownames(geno_df) <- fam$sample.ID      


geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/HSD17B2SNPs.xlsx")

### SRD5A1

library(data.table)
library(SNPassoc)  
library(qqman)     
library(ggplot2)   
library(dplyr)
library(bigsnpr)


snp_readBed("GWAS/SRD5A1.bed")

snp_obj <- snp_attach("GWAS/SRD5A1.rds")

map <- snp_obj$map
head(map)  

fam <- snp_obj$fam
head(fam) 

G <- snp_obj$genotypes

cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)


cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF

allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)
maf


low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID  
all_snps

library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])          
colnames(geno_df) <- map$marker.ID     
rownames(geno_df) <- fam$sample.ID     

geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/SRD5A1SNPs.xlsx")

### SRD5A2

library(data.table)
library(SNPassoc)  
library(qqman)    
library(ggplot2)  
library(dplyr)
library(bigsnpr)


snp_readBed("GWAS/SRD5A2.bed")

snp_obj <- snp_attach("GWAS/SRD5A2.rds")

map <- snp_obj$map
head(map)  


fam <- snp_obj$fam
head(fam) 

G <- snp_obj$genotypes

cat("Performing sample QC...\n")
sample_miss <- rowMeans(is.na(G[]))
summary(sample_miss)

# SNP QC
cat("Performing SNP QC...\n")

snp_miss <- colMeans(is.na(G[]))
keep_snps1 <- which(snp_miss <= 0.05)
summary(snp_miss)
summary(keep_snps1)

# MAF
allele_freq <- colMeans(G[], na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)
maf

low_maf_snps <- names(maf)[maf < 0.01]
low_maf_snps

map <- snp_obj$map

all_snps <- map$marker.ID  
all_snps


library(writexl)
library(openxlsx)

geno_df <- as.data.frame(G[])            
colnames(geno_df) <- map$marker.ID      
rownames(geno_df) <- fam$sample.ID     


geno_df$sample.ID <- rownames(geno_df)

geno_df <- geno_df[, c(ncol(geno_df), 1:(ncol(geno_df)-1))]
write.xlsx(geno_df, file = "GWAS/SRD5A2SNPs.xlsx")


library(openxlsx)

#  Loop through the following genes, load their genotype data

# ESR1 read SNP dataset
geno_df <- read.xlsx("GWAS/SNPs.xlsx")

# ESR2
geno_df <- read.xlsx("GWAS/ESR2SNPs.xlsx")

# AR
geno_df <- read.xlsx("GWAS/ARSNPs.xlsx")

# AKR1C1
geno_df <- read.xlsx("GWAS/AKR1C1SNPs.xlsx")

# CYP19A1
geno_df <- read.xlsx("GWAS/CYP19A1SNPs.xlsx")

# HSD17B2
geno_df <- read.xlsx("GWAS/HSD17B2SNPs.xlsx")

# SRD5A1
geno_df <- read.xlsx("GWAS/SRD5A1SNPs.xlsx")

# SRD5A2
geno_df <- read.xlsx("GWAS/SRD5A2SNPs.xlsx")

# extract incidence bladder cancer dataset from last project
# define as "incident_bladder"

colnames(geno_df)[1] <- "sample.ID"

colnames(incident_bladder)[1] <- "sample.ID"

setdiff(incident_bladder$sample.ID, geno_df$sample.ID)

incident_bladder$sample.ID <- as.character(incident_bladder$sample.ID)
geno_df$sample.ID <- as.character(geno_df$sample.ID)

library(dplyr)
merged_df <- inner_join(incident_bladder, geno_df, by = "sample.ID")

## extract PCs V26-V65, define as "final_snps_PCs"

library(dplyr)

merged_df$sample.ID <- as.character(merged_df$sample.ID)
final_snps_PCs$sample.ID <- as.character(final_snps_PCs$sample.ID)

final_merged <- inner_join(merged_df, final_snps_PCs, by = "sample.ID")

merged_df<-final_merged

### HRT use
merged_df <- merged_df %>%
  mutate(
    HRT_use_full = Ever.used.hormone.replacement.therapy..HRT....Instance.0,
    HRT_use_full = factor(HRT_use_full, 
                          levels = c("No", "Yes", "Do not know", ""), 
                          labels = c("Never", "Ever", "Do not know", "Unknown"))
  )


HRT_count <- table(merged_df$HRT_use_full, useNA = "ifany")
HRT_percent <- prop.table(HRT_count) * 100

data.frame(
  HRT_use = names(HRT_count),
  Count = as.numeric(HRT_count),
  Percent = round(as.numeric(HRT_percent), 2)
)


# filter Do not know and Unknown, 672 subjects in total
merged_df <- merged_df %>%
  filter(!(HRT_use_full %in% c("Do not know", "Unknown")))

table(merged_df$HRT_use_full)

HRT_count <- table(merged_df$HRT_use_full)
HRT_percent <- prop.table(HRT_count) * 100

data.frame(
  HRT_use = names(HRT_count),
  Count = as.numeric(HRT_count),
  Percent = round(as.numeric(HRT_percent), 2)
)

# Rename long variable names to shorter, more readable names
merged_df <- merged_df %>%
  rename(
    assessment_date = Date.of.attending.assessment.centre...Instance.0.x,
    cancer_diagnosis_date = Date.of.cancer.diagnosis...Instance.0,
    age_recruitment = Age.at.recruitment.x)

# Calculate date of birth
merged_df$assessment_date <- as.Date(merged_df$assessment_date)
merged_df$cancer_diagnosis_date <- as.Date(merged_df$cancer_diagnosis_date)

merged_df$DOB <- merged_df$assessment_date -
  merged_df$age_recruitment * 365.25

# Calculate age at cancer diagnosis
merged_df$age_at_diagnosis <- as.numeric(
  (merged_df$cancer_diagnosis_date - merged_df$DOB) / 365.25)

shapiro.test(merged_df$age_at_diagnosis)

# HRT

library(dplyr)

merged_df <- merged_df %>%
  mutate(
    HRT_use = factor(
      HRT_use_full,
      levels = c("Never", "Ever")
    )
  )


merged_df %>%
  group_by(HRT_use) %>%
  summarise(
    N = sum(!is.na(age_at_diagnosis)),
    Median = median(age_at_diagnosis, na.rm = TRUE),
    Q1 = quantile(age_at_diagnosis, 0.25, na.rm = TRUE),
    Q3 = quantile(age_at_diagnosis, 0.75, na.rm = TRUE)
  )


overall_age <- merged_df %>%
  summarise(
    N = sum(!is.na(age_at_diagnosis)),
    Median = median(age_at_diagnosis, na.rm = TRUE),
    Q1 = quantile(age_at_diagnosis, 0.25, na.rm = TRUE),
    Q3 = quantile(age_at_diagnosis, 0.75, na.rm = TRUE)
  )
overall_age

wilcox.test(age_recruitment ~ HRT_use, data = merged_df)

#BMI
library(dplyr)

merged_df <- merged_df %>%
  mutate(
    BMI_category = case_when(
      is.na(Body.mass.index) | Body.mass.index == "" ~ "Unknown",
      Body.mass.index< 25 ~ "Healthy (≤24.9)",
      Body.mass.index>= 25 & Body.mass.index< 30 ~ "Overweight (25–29.9)",
      Body.mass.index>= 30 ~ "Obese (≥30)",
      TRUE ~ "Unknown"
    )
  )

bmi_summary <- merged_df %>%
  group_by(HRT_use, BMI_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HRT_use) %>%
  mutate(prop = n / sum(n) * 100)

chi_result <- chisq.test(table(merged_df$HRT_use, merged_df$BMI_category))

bmi_summary
chi_result$p.value


library(dplyr)

bmi_overall <- merged_df %>%
  group_by(BMI_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  arrange(BMI_category)

bmi_overall

merged_df_bmi <- merged_df %>%
  filter(BMI_category != "Unknown")

chi_result <- chisq.test(table(merged_df_bmi$HRT_use, merged_df_bmi$BMI_category))

chi_result$p.value

### ethnicity

merged_df <- merged_df %>%
  mutate(
    ethnicity_group = case_when(
      Ethnic.background %in% c("British", "Any other white background", "Irish", "White") ~ "White",
      Ethnic.background %in% c("Caribbean", "African", "Any other Black background") ~ "Black",
      Ethnic.background %in% c("Indian", "Any other Asian background", "Chinese", "Bangladeshi", "Pakistani") ~ "Asian",
      Ethnic.background %in% c(
        "Any other mixed background", "White and Asian", "White and Black African", 
        "Asian or Asian British", "Black or Black British", "White and Black Caribbean", "Mixed"
      ) ~ "Mixed",
      Ethnic.background %in% c("Prefer not to answer", "Do not know", "", NA) ~ "Unknown",
      Ethnic.background == "Other ethnic group" ~ "Others",
      TRUE ~ Ethnic.background
    )
  )

ethnicity_summary <- merged_df %>%
  group_by(HRT_use, ethnicity_group) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  arrange(HRT_use, desc(n))

ethnicity_overall <- merged_df %>%
  group_by(ethnicity_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = round(100 * n / sum(n), 2),
         HRT_use = "Overall") %>%
  select(HRT_use, everything())

ethnicity_all <- bind_rows(ethnicity_summary, ethnicity_overall) %>%
  arrange(HRT_use, desc(n))


print(ethnicity_all, n = Inf)

merged_df_ethnicity <- merged_df %>%
  filter(ethnicity_group != "Unknown")

ethnicity_table <- table(merged_df_ethnicity$HRT_use, merged_df_ethnicity$ethnicity_group)

chi_result <- chisq.test(ethnicity_table)

chi_result$p.value

### education

library(dplyr)
library(stringr)
library(purrr)

merged_df <- merged_df %>%
  mutate(education_group = case_when(
    str_detect(Qualifications, regex("None of the above|Prefer not to answer|^$", ignore_case = TRUE)) ~ "Unknown",
    str_detect(Qualifications, regex("CSEs or equivalent", ignore_case = TRUE)) ~ "Less than high school",
    str_detect(Qualifications, regex("O levels/GCSEs or equivalent|A levels/AS levels or equivalent|Other professional qualifications", ignore_case = TRUE)) ~ "High school",
    str_detect(Qualifications, regex("College or University degree|NVQ or HND or HNC or equivalent", ignore_case = TRUE)) ~ "More than high school",
    TRUE ~ "Unknown"
  ))

education_table <- merged_df %>%
  group_by(HRT_use, education_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HRT_use) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  arrange(HRT_use, desc(n))

education_table

education_overall <- merged_df %>%
  group_by(education_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = n / sum(n) * 100,
         HRT_use = "Overall") %>%
  select(HRT_use, everything()) %>%
  arrange(desc(n))

education_overall

library(tidyr)

education_chi <- merged_df %>%
  filter(education_group != "Unknown") %>%
  count(HRT_use, education_group) %>%
  pivot_wider(names_from = education_group, values_from = n, values_fill = 0)

tab <- as.data.frame(education_chi)
rownames(tab) <- tab$HRT_use
tab$HRT_use <- NULL

chisq.test(as.matrix(tab))

### income

library(dplyr)
library(purrr)
library(tidyr)

merged_df <- merged_df %>%
  mutate(
    income_clean = case_when(
      Average.total.household.income.before.tax %in% c("", "Do not know", "Prefer not to answer") ~ "Unknown",
      TRUE ~ Average.total.household.income.before.tax
    )
  )

income_table <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(HRT_use, income_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HRT_use) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  arrange(HRT_use, desc(n))

income_table

income_overall <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(income_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = round(100 * n / sum(n), 2),
         HRT_use = "Overall") %>%
  select(HRT_use, everything()) %>%
  arrange(desc(n))

income_overall

income_chi <- merged_df %>%
  filter(!is.na(HRT_use), income_clean != "Unknown") %>%
  {chisq.test(table(.$HRT_use, .$income_clean))}

income_chi

### smoking

library(dplyr)

merged_df <- merged_df %>%
  mutate(
    Smoking_status_group = ifelse(
      is.na(Smoking.status) | Smoking.status == "",
      "Unknown",
      Smoking.status
    )
  )

smoking_table <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(HRT_use, Smoking_status_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HRT_use) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  arrange(HRT_use, desc(n))

smoking_table

smoking_overall <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(Smoking_status_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    percent = round(100 * n / sum(n), 2),
    HRT_use = "Overall"
  ) %>%
  select(HRT_use, everything()) %>%
  arrange(desc(n))

smoking_overall

smoking_chi <- merged_df %>%
  filter(!is.na(HRT_use), Smoking_status_group != "Unknown") %>%
  {chisq.test(table(.$HRT_use, .$Smoking_status_group))}

smoking_chi

smoking_chi$expected      
smoking_chi$observed   

contribution <- (smoking_chi$observed - smoking_chi$expected)^2 / smoking_chi$expected
contribution


###ocp

library(dplyr)
library(purrr)
library(tidyr)

merged_df <- merged_df %>%
  mutate(
    ocp_use_clean = case_when(
      Ever.taken.oral.contraceptive.pill %in% c("Do not know", "Prefer not to answer", "") ~ "Unknown",
      Ever.taken.oral.contraceptive.pill == "Yes" ~ "Yes",
      Ever.taken.oral.contraceptive.pill == "No" ~ "No",
      TRUE ~ "Unknown"
    )
  )

ocp_table <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(HRT_use, ocp_use_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HRT_use) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(HRT_use, desc(n))

ocp_table

ocp_overall <- merged_df %>%
  filter(!is.na(HRT_use)) %>%
  group_by(ocp_use_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    percent = round(100 * n / sum(n), 1),
    HRT_use = "Overall"
  ) %>%
  select(HRT_use, everything()) %>%
  arrange(desc(n))

ocp_overall

ocp_chi <- merged_df %>%
  filter(!is.na(HRT_use), ocp_use_clean != "Unknown") %>%
  {chisq.test(table(.$HRT_use, .$ocp_use_clean))}

ocp_chi

ocp_contribution <- (ocp_chi$observed - ocp_chi$expected)^2 / ocp_chi$expected
ocp_contribution

### ocp duration

library(dplyr)

merged_df <- merged_df %>%
  mutate(
    ocp_age_start = as.numeric(Age.started.oral.contraceptive.pill),
    ocp_age_end = as.numeric(Age.when.last.used.oral.contraceptive.pill),
    ocp_duration = case_when(
      Age.started.oral.contraceptive.pill %in% c("Do not know", "Prefer not to answer", "Still taking the pill", "") ~ NA_real_,
      Age.when.last.used.oral.contraceptive.pill %in% c("Do not know", "Prefer not to answer", "Still taking the pill", "") ~ NA_real_,
      TRUE ~ ocp_age_end - ocp_age_start
    )
  )

ocp_duration_summary <- merged_df %>%
  filter(!is.na(ocp_duration), !is.na(HRT_use)) %>%
  group_by(HRT_use) %>%
  summarise(
    n = n(),
    median_duration = median(ocp_duration, na.rm = TRUE),
    IQR_lower = quantile(ocp_duration, 0.25, na.rm = TRUE),
    IQR_upper = quantile(ocp_duration, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(ocp_duration_summary)

ocp_duration_summary_overall <- merged_df %>%
  filter(!is.na(ocp_duration)) %>%
  summarise(
    n = n(),
    median_duration = median(ocp_duration, na.rm = TRUE),
    IQR_lower = quantile(ocp_duration, 0.25, na.rm = TRUE),
    IQR_upper = quantile(ocp_duration, 0.75, na.rm = TRUE)
  )

print(ocp_duration_summary_overall)

ocp_duration_wilcox <- merged_df %>%
  filter(!is.na(ocp_duration), !is.na(HRT_use)) 

if(length(unique(ocp_duration_wilcox$HRT_use)) > 1){
  wilcox_test <- wilcox.test(ocp_duration ~ HRT_use, data = ocp_duration_wilcox)
  wilcox_test$p.value
} else {
  NA
}


library(survival)
### collinearity check

cox_data <- merged_df %>%
  filter(
    Ever.taken.oral.contraceptive.pill != "Prefer not to answer",
    Smoking.status != "Prefer not to answer"
  )

library(car)

covariates_df <- cox_data %>%
  select(HRT_use, age_at_diagnosis, 
         Ever.taken.oral.contraceptive.pill, 
         Smoking.status) %>%
  mutate(across(where(is.character), as.factor))

lm_fit <- lm(as.numeric(HRT_use) ~ age_at_diagnosis + 
               Ever.taken.oral.contraceptive.pill+ 
               Smoking.status, data = covariates_df)

vif(lm_fit)

library(survival)
library(dplyr)

### Overall Survival (OS)

### OS

# Ensure date formats are correct
merged_df$diagnosis_date <- as.Date(merged_df$Date.of.cancer.diagnosis)
merged_df$death_date <- as.Date(merged_df$Date.of.death)

# Follow-up end date
followup_end <- as.Date("2024-12-31")

# Select end date: death date if available, otherwise follow-up end date
merged_df$end_date <- merged_df$death_date
merged_df$end_date[is.na(merged_df$end_date)] <- followup_end

# Calculate overall survival time (years)
merged_df$OS_time <- as.numeric(
  merged_df$end_date - merged_df$diagnosis_date
) / 365

# Check
summary(merged_df$OS_time)

merged_df$OS_event <- ifelse(!is.na(merged_df$death_date), 1, 0)

table(merged_df$OS_event)

surv_obj <- Surv(time = merged_df$OS_time, event = merged_df$OS_event)

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# HRT KM

table(merged_df$HRT_use, useNA = "ifany")


library(survival)
library(survminer)

km_fit <- survfit(surv_obj ~ HRT_use, data = merged_df)

ggsurvplot(
  km_fit,
  data = merged_df,
  risk.table = TRUE,          
  pval = TRUE,                
  conf.int = TRUE,            
  xlab = "Time (years)",
  ylab = "Overall Survival Probability",
  legend.title = "HRT Use",
  legend.labs = c("Never", "Ever"),
  palette = c("#E7B800", "#2E9FDF")
)

library(dplyr)

hrt_summary <- merged_df %>%
  group_by(HRT_use) %>%   # HRT_use = Yes / No
  summarise(
    N = n(),
    events = sum(OS_event, na.rm = TRUE)
  )

hrt_summary


# 5 years
surv_5yr <- summary(km_fit, times = 5)

# output
surv_5yr_table <- data.frame(
  HRT_use = surv_5yr$strata,
  time = surv_5yr$time,
  survival = round(surv_5yr$surv, 3),
  lower_CI = round(surv_5yr$lower, 3),
  upper_CI = round(surv_5yr$upper, 3)
)

surv_5yr_table


surv_obj <- Surv(time = merged_df$OS_time, event = merged_df$OS_event)

# median + 95%CI
fit_overall <- survfit(surv_obj ~ 1)

summary(fit_overall)$table[c("median", "0.95LCL", "0.95UCL")]

### HRT cox regression


# remove OCP and Smoking's "Prefer not to answer"
cox_data <- merged_df %>%
  filter(
    Ever.taken.oral.contraceptive.pill!= "Prefer not to answer",
    Smoking.status!= "Prefer not to answer"
  )

cox_fit <- coxph(
  Surv(OS_time, OS_event) ~ 
    HRT_use + 
    age_at_diagnosis + 
    Ever.taken.oral.contraceptive.pill + 
    Smoking.status,
  data = cox_data
)

summary(cox_fit)


ph_test <- cox.zph(cox_fit)

ph_test

####HRT duration
library(dplyr)
library(survival)

merged_df$Age.started.hormone.replacement.therapy[merged_df$Age.started.hormone.replacement.therapy== ""] <- NA
merged_df$Age.last.used.hormone.replacement.therapy[merged_df$Age.last.used.hormone.replacement.therapy== ""] <- NA

merged_df$Age.started.hormone.replacement.therapy<- as.numeric(merged_df$Age.started.hormone.replacement.therapy)
merged_df$Age.last.used.hormone.replacement.therapy<- as.numeric(merged_df$Age.last.used.hormone.replacement.therapy)

merged_df$HRT_duration <- merged_df$Age.last.used.hormone.replacement.therapy -
  merged_df$Age.started.hormone.replacement.therapy

summary(merged_df$HRT_duration)

num_missing <- sum(is.na(merged_df$HRT_duration))
percent_missing <- mean(is.na(merged_df$HRT_duration)) * 100

cat("HRT duration missing:", num_missing, "subjects (", round(percent_missing, 2), "%)\n")

covariates_df <- merged_df %>%
  select(HRT_duration, age_at_diagnosis, 
         Ever.taken.oral.contraceptive.pill, 
         Smoking.status) %>%
  mutate(across(where(is.character), as.factor))

lm_fit <- lm(HRT_duration ~ age_at_diagnosis + 
               Ever.taken.oral.contraceptive.pill + 
               Smoking.status, data = covariates_df)

vif(lm_fit)

surv_obj <- Surv(time = merged_df$OS_time, event = merged_df$OS_event)



# remove  OCP and Smoking 's "Prefer not to answer"
cox_data <- merged_df %>%
  filter(
    Ever.taken.oral.contraceptive.pill != "Prefer not to answer",
    Smoking.status != "Prefer not to answer"
  )

cox_fit <- coxph(
  Surv(OS_time, OS_event) ~ 
    HRT_duration + 
    age_at_diagnosis + 
    Ever.taken.oral.contraceptive.pill+ 
    Smoking.status,
  data = cox_data
)

summary(cox_fit)

ph_test <- cox.zph(cox_fit)

ph_test

# Loop through the following genes, load their genotype data

#ESR1

snps <- c("rs1293940", "rs851995", "rs2941740", "rs73005768", "rs12526447",
          "rs7772579", "rs1415194", "rs1999805", "rs2152750", "rs2504069",
          "rs2504063", "rs17081685", "rs538098", "rs139960913", "rs2077647",
          "rs9371557", "rs9340789", "rs11155813", "rs6902771", "rs4986936",
          "rs2234693", "rs9340799", "rs78796841", "rs2175898", "rs35882398",
          "rs9322341", "rs12204714", "rs9383951", "rs12212176", "rs722208",
          "rs926779", "rs9341086", "rs2295192")

#ESR2
snps <- c("rs10137972", "rs2781377", "rs12881815", "rs12883524", "rs36215895",
          "rs1152591", "rs1152590", "rs12435857", "rs117171778")

# AR
snps <- c("rs5031002")

# AKR1C1
snps <- c("rs2904804", "rs77045180", "rs76207412")

# CYP19A1
snps <- c("rs4646", "rs700519", "rs2899472", "rs17523014", "rs12439137",
          "rs4775935", "rs2414095", "rs700518", "rs16964211", "rs10459592",
          "rs2305707", "rs3751591", "rs7174997")

# HSD17B2
snps <- c("rs72829163", "rs200540657", "rs79488654")

# SRD5A1
snps <- c("rs472402")

# SRD5A2
snps <- c("rs522638", "rs523349")

# PCs
pcs <- paste0("V", 26:35)

results <- data.frame(
  SNP = snps,
  HR = NA,
  lower95 = NA,
  upper95 = NA,
  p.value = NA
)

for (i in seq_along(snps)) {
  snp_name <- snps[i]
  
  covariates <- paste(c("age_at_diagnosis", pcs), collapse = " + ")
  fml <- as.formula(paste0("surv_obj ~ ", snp_name, " + ", covariates))
  
  fit <- coxph(fml, data = merged_df)
  sum_fit <- summary(fit)
  
  results$HR[i] <- sum_fit$coef[1, "exp(coef)"]
  results$lower95[i] <- sum_fit$conf.int[1, "lower .95"]
  results$upper95[i] <- sum_fit$conf.int[1, "upper .95"]
  results$p.value[i] <- sum_fit$coef[1, "Pr(>|z|)"]
}

results

##PH

ph_results <- data.frame(
  SNP = snps,
  GLOBAL_p = NA
)

for (i in seq_along(snps)) {
  snp_name <- snps[i]
  
  covariates <- paste(c("age_at_diagnosis", pcs), collapse = " + ")
  fml <- as.formula(paste0("surv_obj ~ ", snp_name, " + ", covariates))
  
  fit <- coxph(fml, data = merged_df)
  
  ph_test <- cox.zph(fit)
  
  ph_results$GLOBAL_p[i] <- ph_test$table["GLOBAL", "p"]
}

ph_results

###GRS

##GRSOS

library(openxlsx)
library(dplyr)

ESR1<- read.xlsx("GWAS/SNPs.xlsx")
CYP19A1<- read.xlsx("GWAS/CYP19A1SNPs.xlsx")

ESR1_sel   <- ESR1   %>% select(sample.ID, rs2295192)
CYP19A1_sel <- CYP19A1  %>% select(sample.ID, rs700519,rs16964211)

snps_merged <- ESR1_sel %>%
  inner_join(CYP19A1_sel, by = "sample.ID")

head(snps_merged)


beta_rs2295192 <- log(2.4412454)
beta_rs700519 <- log(1.8478048)
beta_rs16964211 <- log(1.5298841)

beta_rs2295192
beta_rs700519
beta_rs16964211

merged_full <- inner_join(merged_df, snps_merged, by = "sample.ID")

# weighted
merged_full$GRS_weighted <- 
  beta_rs2295192 * merged_full$rs2295192 +
  beta_rs700519 * merged_full$rs700519 +
  beta_rs16964211 * merged_full$rs16964211

# standardized（mean=0，SD=1）
merged_full$GRS_std <- as.numeric(scale(merged_full$GRS_weighted))

# Shapiro-Wilk
shapiro.test(merged_full$GRS_std)   

baseline_grs <- merged_full %>%
  group_by(HRT_use) %>%
  summarise(
    n = n(),
    median_PRS = median(GRS_std, na.rm = TRUE),
    Q1 = quantile(GRS_std, 0.25, na.rm = TRUE),
    Q3 = quantile(GRS_std, 0.75, na.rm = TRUE)
  )
baseline_grs

HRT_grs_summary <- merged_full %>%
  group_by(HRT_use) %>%
  summarise(
    n = n(),
    mean_GRS = mean(GRS_std, na.rm = TRUE),
    sd_GRS = sd(GRS_std, na.rm = TRUE)
  )

HRT_grs_summary

wilcox.test(GRS_std ~ HRT_use,
            data = merged_full %>% filter(HRT_use %in% c("Never", "Ever")))
overall_GRS_summary <- merged_full %>%
  summarise(
    n = sum(!is.na(GRS_std)),          
    mean_GRS = mean(GRS_std, na.rm = TRUE),
    sd_GRS = sd(GRS_std, na.rm = TRUE)
  )

overall_GRS_summary

library(survival)

surv_os <- Surv(time = merged_full$OS_time, event = merged_full$OS_event)

pcs <- paste0("V", 26:35) 
formula_str <- paste("surv_os ~ GRS_std + age_at_diagnosis +", paste(pcs, collapse = " + "))
cox_model <- coxph(as.formula(formula_str), data = merged_full)

summary(cox_model)

ph_test <- cox.zph(cox_model)

ph_test

plot(ph_test)

### Interaction

####GRS*

library(survival)
library(dplyr)
library(tidyr)

merged_full <- merged_full %>%
  filter(HRT_use %in% c("Never", "Ever")) %>%
  mutate(HRT_use = factor(HRT_use, levels = c("Never", "Ever")))

merged_full_cox_OS <- merged_full %>%
  select(OS_time, OS_event, GRS_std, age_at_diagnosis, HRT_use, all_of(pcs)) %>%
  drop_na()

surv_os <- Surv(time = merged_full_cox_OS$OS_time,
                event = merged_full_cox_OS$OS_event)

formula_os <- paste(
  "surv_os ~ GRS_std * HRT_use + age_at_diagnosis +",
  paste(pcs, collapse = " + ")
)

cox_interaction_OS <- coxph(as.formula(formula_os),
                            data = merged_full_cox_OS)

summary(cox_interaction_OS)


ph_test <- cox.zph(cox_interaction_OS)

ph_test


### CSS

merged_df$CS_event <- ifelse(grepl("^C67", merged_df$Cause.of.deathICD.10), 1, 0)

followup_end <- as.Date("2024-12-31")

merged_df$CS_event <- ifelse(!is.na(merged_df$Date.of.death) & 
                               grepl("^C67", merged_df$Cause.of.deathICD.10),
                             1, 0)


merged_df$CS_end_date <- ifelse(merged_df$CS_event == 1,
                                merged_df$Date.of.death,
                                followup_end)

merged_df$CS_end_date <- as.character(merged_df$CS_end_date)

merged_df$CS_end_date[merged_df$CS_end_date == "20088"] <- NA

merged_df$CS_end_date <- as.Date(merged_df$CS_end_date, format="%Y/%m/%d")

merged_df$CS_end_date[is.na(merged_df$CS_end_date)] <- followup_end

merged_df$CS_end_date <- as.Date(merged_df$CS_end_date)
merged_df$diagnosis_date <- as.Date(merged_df$Date.of.cancer.diagnosis)

merged_df$CS_time <- as.numeric(merged_df$CS_end_date - merged_df$diagnosis_date) / 365

summary(merged_df$CS_time)
table(merged_df$CS_event)

surv_cs <- Surv(time = merged_df$CS_time, event = merged_df$CS_event)

# HRT KM CSS
km_fit_css <- survfit(surv_cs ~ HRT_use, data = merged_df)


ggsurvplot(
  km_fit_css,
  data = merged_df,
  risk.table = TRUE,         
  pval = TRUE,                
  conf.int = TRUE,       
  xlab = "Time (years)",
  ylab = "Cancer-Specific Survival Probability",
  legend.title = "HRT Use",
  legend.labs = c("Never", "Ever"),
  palette = c("#E7B800", "#2E9FDF")
)

hrt_summary_css <- merged_df %>%
  group_by(HRT_use) %>%
  summarise(
    N = n(),
    events = sum(CS_event, na.rm = TRUE)
  )

hrt_summary_css

median_surv <- summary(km_fit_css)$table[, c("median", "0.95LCL", "0.95UCL")]
median_surv <- as.data.frame(median_surv)
median_surv$HRT_use <- rownames(median_surv)
rownames(median_surv) <- NULL
median_surv <- median_surv %>% select(HRT_use, everything())
median_surv

surv_5yr <- summary(km_fit_css, times = 5)$surv
surv_5yr_lower <- summary(km_fit_css, times = 5)$lower
surv_5yr_upper <- summary(km_fit_css, times = 5)$upper
hr_labels <- names(km_fit_css$strata)

surv_5yr_df <- data.frame(
  HRT_use = gsub("HRT_use=", "", hr_labels),
  surv_5yr = surv_5yr,
  lower_95CI = surv_5yr_lower,
  upper_95CI = surv_5yr_upper
)
surv_5yr_df

km_summary <- summary(km_fit_css)

km_summary

# cox regression
library(survival)
library(dplyr)

merged_df_clean <- merged_df %>%
  filter(
    Ever.taken.oral.contraceptive.pill!= "Prefer not to answer",
    Smoking.status!= "Prefer not to answer"
  )

# factor
merged_df_clean$Ever.taken.oral.contraceptive.pill<- 
  factor(merged_df_clean$Ever.taken.oral.contraceptive.pill)

merged_df_clean$Smoking.status<- 
  factor(merged_df_clean$Smoking.status)

cox_fit_css <- coxph(
  surv_cs ~ 
    HRT_use + 
    age_at_diagnosis + 
    Ever.taken.oral.contraceptive.pill + 
    Smoking.status,
  data = merged_df_clean
)

summary(cox_fit_css)

ph_test <- cox.zph(cox_fit_css)

ph_test


### HRT duration
cox_fit_css <- coxph(
  surv_cs ~ 
    HRT_duration + 
    age_at_diagnosis + 
    Ever.taken.oral.contraceptive.pill + 
    Smoking.status,
  data = merged_df
)

summary(cox_fit_css)

ph_test <- cox.zph(cox_fit_css)

ph_test

#Loop through the following genes, load their genotype data
#ESR1

snps <- c("rs1293940", "rs851995", "rs2941740", "rs73005768", "rs12526447",
          "rs7772579", "rs1415194", "rs1999805", "rs2152750", "rs2504069",
          "rs2504063", "rs17081685", "rs538098", "rs139960913", "rs2077647",
          "rs9371557", "rs9340789", "rs11155813", "rs6902771", "rs4986936",
          "rs2234693", "rs9340799", "rs78796841", "rs2175898", "rs35882398",
          "rs9322341", "rs12204714", "rs9383951", "rs12212176", "rs722208",
          "rs926779", "rs9341086", "rs2295192")

#ESR2
snps <- c("rs10137972", "rs2781377", "rs12881815", "rs12883524", "rs36215895",
          "rs1152591", "rs1152590", "rs12435857", "rs117171778")

# AR
snps <- c("rs5031002")

# AKR1C1
snps <- c("rs2904804", "rs77045180", "rs76207412")

# CYP19A1
snps <- c("rs4646", "rs700519", "rs2899472", "rs17523014", "rs12439137",
          "rs4775935", "rs2414095", "rs700518", "rs16964211", "rs10459592",
          "rs2305707", "rs3751591", "rs7174997")

# HSD17B2
snps <- c("rs72829163", "rs200540657", "rs79488654")

# SRD5A1
snps <- c("rs472402")

# SRD5A2
snps <- c("rs522638", "rs523349")

# PCs 
pcs <- paste0("V", 26:35)

results <- data.frame(
  SNP = snps,
  HR = NA,
  lower95 = NA,
  upper95 = NA,
  p.value = NA
)

# analyses CSS
for (i in seq_along(snps)) {
  snp_name <- snps[i]
  
  covariates <- paste(c("age_at_diagnosis", pcs), collapse = " + ")
  fmls <- as.formula(paste0("surv_cs ~ ", snp_name, " + ", covariates))
  
  fit <- coxph(fmls, data = merged_df)
  sum_fit <- summary(fit)
  
  results$HR[i] <- sum_fit$coef[1, "exp(coef)"]
  results$lower95[i] <- sum_fit$conf.int[1, "lower .95"]
  results$upper95[i] <- sum_fit$conf.int[1, "upper .95"]
  results$p.value[i] <- sum_fit$coef[1, "Pr(>|z|)"]
}

results

##PH

ph_results_css <- data.frame(
  SNP = snps,
  GLOBAL_p = NA
)

for (i in seq_along(snps)) {
  snp_name <- snps[i]
  
  covariates <- paste(c("age_at_diagnosis", pcs), collapse = " + ")
  fml <- as.formula(paste0("surv_cs ~ ", snp_name, " + ", covariates))
  
  fit <- coxph(fml, data = merged_df)
  
  ph_test <- cox.zph(fit)
  
  ph_results_css$GLOBAL_p[i] <- ph_test$table["GLOBAL", "p"]
}

ph_results_css

###GRS

##GRSCSS

library(openxlsx)
library(dplyr)


ESR2   <- read.xlsx("GWAS/ESR2SNPs.xlsx")
SRD5A2    <- read.xlsx("GWAS/SRD5A2SNPs.xlsx")

ESR2_sel   <- ESR2   %>% select(sample.ID, rs12881815)
SRD5A2_sel <- SRD5A2    %>% select(sample.ID, rs522638,rs523349)


snps_merged <- ESR2_sel %>%
  inner_join(SRD5A2_sel, by = "sample.ID")

head(snps_merged)


beta_rs12881815 <- log(2.145591e-01)
beta_rs522638 <- log(0.5992727)
beta_rs523349 <- log(0.5655452)

beta_rs12881815
beta_rs522638
beta_rs523349

merged_full <- inner_join(merged_df, snps_merged, by = "sample.ID")

merged_full$GRS_weighted <- 
  beta_rs12881815 * merged_full$rs12881815 +
  beta_rs522638 * merged_full$rs522638 +
  beta_rs523349 * merged_full$rs523349


merged_full$GRS_std <- as.numeric(scale(merged_full$GRS_weighted))

# GRS
shapiro.test(merged_full$GRS_std)  


baseline_grs <- merged_full %>%
  group_by(HRT_use) %>%
  summarise(
    n = n(),
    median_PRS = median(GRS_std, na.rm = TRUE),
    Q1 = quantile(GRS_std, 0.25, na.rm = TRUE),
    Q3 = quantile(GRS_std, 0.75, na.rm = TRUE)
  )
baseline_grs


HRT_grs_summary <- merged_full %>%
  group_by(HRT_use) %>%
  summarise(
    n = n(),
    mean_GRS = mean(GRS_std, na.rm = TRUE),
    sd_GRS = sd(GRS_std, na.rm = TRUE)
  )

HRT_grs_summary

wilcox.test(GRS_std ~ HRT_use,
            data = merged_full %>% filter(HRT_use %in% c("Never", "Ever")))

overall_GRS_summary <- merged_full %>%
  summarise(
    n = sum(!is.na(GRS_std)),         
    mean_GRS = mean(GRS_std, na.rm = TRUE),
    sd_GRS = sd(GRS_std, na.rm = TRUE)
  )

overall_GRS_summary


library(survival)

merged_full$CS_event <- ifelse(grepl("^C67", merged_full$Cause.of.deathICD.10), 1, 0)

followup_end <- as.Date("2024-12-31")

merged_full$CS_event <- ifelse(
  !is.na(merged_full$Date.of.death) &
    grepl("^C67", merged_full$Cause.of.deathICD.10),
  1, 0
)

merged_full$CS_end_date <- ifelse(
  merged_full$CS_event == 1,
  merged_full$Date.of.death,
  followup_end
)


merged_full$CS_end_date <- as.character(merged_full$CS_end_date)
merged_full$CS_end_date[merged_full$CS_end_date == "20088"] <- NA

merged_full$CS_end_date <- as.Date(merged_full$CS_end_date, format="%Y/%m/%d")

merged_full$CS_end_date[is.na(merged_full$CS_end_date)] <- followup_end

merged_full$CS_end_date <- as.Date(merged_full$CS_end_date)
merged_full$diagnosis_date<- as.Date(merged_full$Date.of.cancer.diagnosis)

merged_full$CS_time <- as.numeric(
  merged_full$CS_end_date - merged_full$diagnosis_date
) / 365


summary(merged_full$CS_time)
table(merged_full$CS_event)

library(survival)
surv_cs <- Surv(time = merged_full$CS_time, event = merged_full$CS_event)


pcs <- paste0("V", 26:35)  
formula_str <- paste("surv_cs ~ GRS_std + age_at_diagnosis +", paste(pcs, collapse = " + "))
cox_model <- coxph(as.formula(formula_str), data = merged_full)

summary(cox_model)


ph_test <- cox.zph(cox_model)

ph_test

#  Schoenfeld plot
plot(ph_test)


### Interaction
####GRS

library(survival)
library(dplyr)

merged_full <- merged_full %>%
  filter(HRT_use %in% c("Never", "Ever")) %>%
  mutate(HRT_use = factor(HRT_use, levels = c("Never", "Ever")))

library(tidyr)
merged_full_cox <- merged_full %>%
  select(CS_time, CS_event, GRS_std, age_at_diagnosis, HRT_use, all_of(pcs)) %>%
  drop_na()

surv_cs <- Surv(time = merged_full_cox$CS_time, event = merged_full_cox$CS_event)

formula_str <- paste(
  "surv_cs ~ GRS_std * HRT_use + age_at_diagnosis +",
  paste(pcs, collapse = " + ")
)

cox_interaction <- coxph(as.formula(formula_str), data = merged_full_cox)

summary(cox_interaction)

ph_test <- cox.zph(cox_interaction)

ph_test





