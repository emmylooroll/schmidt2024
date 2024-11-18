install.packages("R.utils")
# Load the required packages (install any that are missing)
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("ieugwasr", quietly = TRUE)) install.packages("ieugwasr")
if (!requireNamespace("coloc", quietly = TRUE)) install.packages("coloc")

library(usethis)
library(openxlsx)
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(plyr)
library(stringr)
library(coloc)


setwd("xx")

usethis::edit_r_environ()

ieugwasr::get_opengwas_jwt()
ieugwasr::user()

gwasinfo <- gwasinfo()

#List of three independent genome-wide significant SNPs in Blauw et al (2018) paper (Table 2)
CETP <-read_exposure_data(
  "exposure file.csv",
  clump = FALSE,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "Effect allele",
  other_allele_col = "Other allele",
  chr_col = "CHR",
  pos_col = "POS",
  pval_col = "P"
)

##manually putting in p va

#Full list of CETP SNPs in Blauw et al (2018) paper - all on CHR16, not pruned
CETP_full <-read_exposure_data(
  "cetp_gws_hits.csv",
  clump = FALSE,
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "PVALUE",
  chr_col = "CHR",
  pos_col = "POS"
)

#Clumping full CETP SNP list
cetp_clumped <- clump_data(
  CETP_full,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR",
  bfile = NULL,
  plink_bin = NULL
)
                
                         
##LOOKING UP the 3 CETP SNPS IN various AD GWAS DATASETS AND RUNNING MR##
##PLEASE NOTE THAT CETP INHIBITORS LOWER CETP, BUT THE MR WILL RUN PER UNIT INCREASE IN CETP.
##RESULTING COEFFICIENTS NEED TO BE FLIPPED TO REPRESENT DRUG EFFECT
#LAMBERT
lambert <-associations(variants=CETP$SNP, id=c("ieu-a-297"))
lambert <- format_data(
  lambert,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR WITH LAMBERT
dat <- harmonise_data(CETP, lambert, action = 2)
lambert_mr_result <- mr(dat)
lambert_mr_result
lambert_mr_result$beta_flipped <- -lambert_mr_result$b  #Flipping log(OR) to reflect inhibiting effect
write.csv(lambert_mr_result, "lambert_results.csv")

#KUNKLE
kunkle <-associations(variants=CETP$SNP, id=c("ieu-b-2"))
kunkle <- format_data(
  kunkle,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR WITH KUNKLE
dat <- harmonise_data(CETP, kunkle, action = 2)
kunkle_mr_result <- mr(dat)
kunkle_mr_result
kunkle_mr_result$beta_flipped <- -kunkle_mr_result$b  #Flipping log(OR) to reflect inhibiting effect
write.csv(kunkle_mr_result, "kunkle_results.csv")

#JANSEN
setwd("xx")
jansen <- read.table("AD_sumstats_Jansenetal_2019sept.txt", header=TRUE)
jansen <- format_data(
  jansen,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "BP"
)

#RUNNING MR WITH JANSEN
dat <- harmonise_data(CETP, jansen, action = 2)
jansen_mr_result <- mr(dat)
jansen_mr_result
jansen_mr_result$beta_flipped <- -jansen_mr_result$b #Flipping log(OR) to reflect inhibiting effect
setwd("xx")
write.csv(jansen_mr_result, "jansen_results.csv")


#BELLENGUEZ
bellenguez <-associations(variants=CETP$SNP, id=c("ebi-a-GCST90027158"))
bellenguez <- format_data(
  bellenguez,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR WITH BELLENGUEZ
dat <- harmonise_data(CETP, bellenguez, action = 2)
bellenguez_mr_result <- mr(dat)
bellenguez_mr_result
bellenguez_mr_result$beta_flipped <- -bellenguez_mr_result$b #Flipping log(OR) to reflect inhibiting effect
write.csv(bellenguez_mr_result, "bellenguez_results.csv")



#SCALING TO HDL, LDL & TG USING WILLER ET AL 2013 SUMSTATS
####PLEASE NOTE THAT CETP INHIBITORS LOWER CETP, RAISE HDL, AND LOWER LDL AND TG. ALL MR WILL RUN PER UNIT INCREASE IN EXPOSURE
##RESULTING COEFFICIENTS FOR CETP, LDL AND TG NEED TO BE FLIPPED TO REPRESENT DRUGS  EFFECT ON DECREASING THESE TRAITS

#HDL
cetp_hdl <-associations(variants=CETP$SNP, id=c("ieu-a-299"))
cetp_hdl <- format_data(
  cetp_hdl,
  type = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR SCALED TO HDL
dat <- harmonise_data(cetp_hdl, bellenguez, action = 2)
hdl_bellenguez_mr_result <- mr(dat)
hdl_bellenguez_mr_result #DO NOT FLIP THIS COEFFICIENT - NEEDS TO REFLECT INCREASE IN HDL
write.csv(hdl_bellenguez_mr_result, "hdl_bellenguez_mr_result.csv")

#LDL
cetp_ldl <-associations(variants=CETP$SNP, id=c("ieu-a-300"))
cetp_ldl <- format_data(
  cetp_ldl,
  type = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR SCALED TO LDL
dat <- harmonise_data(cetp_ldl, bellenguez, action = 2)
ldl_bellenguez_mr_result <- mr(dat)
ldl_bellenguez_mr_result
ldl_bellenguez_mr_result$beta_flipped <- -ldl_bellenguez_mr_result$b  # Flipping log(OR) to reflect lowering of LDL by drug
write.csv(ldl_bellenguez_mr_result, "ldl_bellenguez_results.csv")

#TG
cetp_tg <-associations(variants=CETP$SNP, id=c("ieu-a-302"))
cetp_tg <- format_data(
  cetp_tg,
  type = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos"
)

#RUNNING MR SCALED TO TG
dat <- harmonise_data(cetp_tg, bellenguez, action = 2)
tg_bellenguez_mr_result <- mr(dat)
tg_bellenguez_mr_result
tg_bellenguez_mr_result$beta_flipped <- -tg_bellenguez_mr_result$b  # Flipping log(OR) to reflect lowering of TG by drug
write.csv(tg_bellenguez_mr_result, "tg_bellenguez_results.csv")


###RUNNING MR WITH VAD AS OUTCOME FROM FONGANG ET AL 2024
setwd("xx")
vad <- read.table("MEGAVCID_VaD_euro.gz", fill=TRUE, header=TRUE)
vad <- format_data(
  vad,
  type = "outcome",
  snp_col = "rsID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P.value",
  chr_col = "CHR",
  pos_col = "POS"
)

dat <- harmonise_data(CETP, vad, action = 2)
vad_mr_result <- mr(dat)
vad_mr_result
vad_mr_result$beta_flipped <- -vad_mr_result$b  # Flipping log(OR)
setwd("C:/Users/epela/OneDrive - University of Bristol/MIGRATED O DRIVE/Responses to published papers/Schmidt et al - CETP dementia paper/Data")
write.csv(vad_mr_result, "vad_results.csv")

###RUNNING MR WITH ACD AS OUTCOME FROM FONGANG ET AL 2024
setwd("xx")
acd <- read.table("MEGAVCID_ACD_euro.gz", fill=TRUE, header=TRUE)
acd <- format_data(
  acd,
  type = "outcome",
  snp_col = "rsID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P.value",
  chr_col = "CHR",
  pos_col = "POS"
)

dat <- harmonise_data(CETP, acd, action = 2)
acd_mr_result <- mr(dat)
acd_mr_result
acd_mr_result$beta_flipped <- acd_mr_result$b  # Flipping log(OR)
setwd("xx")
write.csv(acd_mr_result, "acd_results.csv")

