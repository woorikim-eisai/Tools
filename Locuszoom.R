#Author: Woori Kim
#Last update: 2024-05-14
#This R script creates a locuszoom plot using locuszoomr package
#Input format: chr, pos (hg19), rsid and p 
#Example
#chrom,pos,rsid,p
#2,113473611,2:113473611_AT_A,0.0001
#2,113473650,rs138922943,0.5463

library(locuszoomr); library(EnsDb.Hsapiens.v75); library(data.table)

# STEP1 ########################################################################
# Prepare input for your analysis

# Example 
Gene="PAX8"
Pop="EUR"
Filename.Input=paste0("/mnt/zfs/P201_UKB_pQTL/results/PAX8/region/plot/input/results_PAX8_region_IL1RN_P18510_OID20700_v1_Inflammation_rsID.txt")
Filename.Output=paste0("/mnt/zfs/P201_UKB_pQTL/results/PAX8/region/plot/output/locuszoom_PAX8_region_IL1RN_P18510_OID20700_v1_Inflammation_EUR.pdf")

# STEP2 ########################################################################
# Create locus object for plotting

# Read input 
dat <- fread(Filename.Input)

# Create locus object for plotting: gene region +/-500kb 
loc <- locus(gene = Gene, data=dat, flank = 5e5, ens_db="EnsDb.Hsapiens.v75")

# Add LD: to do this, get the token in here: https://ldlink.nih.gov/?tab=apiaccess
#it queries LDlink (https://ldlink.nci.nih.gov/) via the LDlinkR package.
#pop: a 1000 Genomes Project population. multiple populations allowed.
loc <- link_LD(loc, 
               pop = Pop,
               r2d = "r2", 
               genome_build="grch37",
               token = "d67110dd7bd9")

# Add recombination rate
#loc <- link_recomb(loc, 
#                   genome = "hg19")
#retrieving UCSC recombination rate takes time and sometimes get error

#use downloaded the whole recombination rate track
#Caution: Used Combined Population
#there are other track for CEU and YRI.
recomb.hg19 <- import.bw("/mnt/zfs/P201_UKB_pQTL/data/hapMapRelease24CombinedRecombMap.bw")
loc <- link_recomb(loc, 
                   recomb = recomb.hg19)

#summary locus object
summary(loc)

# STEP3.1 ########################################################################
# Plot locuszoom plot

pdf(Filename.Output)
locus_plot(loc, 
           labels = c("index"), #can add multiple SNPs
           label_x = c(4, -5), 
           maxrows = 2, 
           filter_gene_biotype = 'protein_coding',
           highlight = Gene)
dev.off()

# Plot locuszoom plotly
locus_plotly(loc,
             maxrows = 2,
             filter_gene_biotype = 'protein_coding')

# STEP3.2 ########################################################################
# Plot multiple locuszoom plot using for loop 
Assay=c("IL1RN_P18510_OID20700_v1_Inflammation", "IL1R1_P14778_OID21116_v1_Neurology", "HEPACAM2_A8MVW5_OID31151_v1_Neurology_II", "CALCA_P01258_OID20983_v1_Neurology", "IL36G_Q9NZH8_OID30551_v1_Inflammation_II", "MMP1_P03956_OID20672_v1_Inflammation", "IFIT3_O14879_OID31283_v1_Oncology_II", "KLK15_Q9H2R5_OID31427_v1_Oncology_II", "TG_P01266_OID31359_v1_Oncology_II", "DPP10_Q8N608_OID20527_v1_Inflammation", "UMOD_P07911_OID20237_v1_Cardiometabolic")

for (i in 1:length(Assay)) {

  # make object for Filename.Input and Filename.Output
  Filename.Input=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/input/results_", Gene, "_region_", Assay[i], "_rsID.txt")
  Filename.Output=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/output/locuszoom_", Gene, "_region_", Assay[i], "_", Pop, ".pdf")
  
  # read input 
  dat <- fread(Filename.Input)
  
  # create locus object for plotting: gene region +/-500kb 
  loc <- locus(gene = Gene, data=dat, flank = 5e5, ens_db="EnsDb.Hsapiens.v75")
  
  #add LD
  loc <- link_LD(loc, 
                 pop = Pop,
                 r2d = "r2", 
                 genome_build="grch37",
                 token = "d67110dd7bd9")
  
  #add recombination rate
  loc <- link_recomb(loc, 
                     genome = "hg19")
  
  #summary locus object
  summary(loc)
  
  #plot locuszoom plot
  pdf(Filename.Output)
  locus_plot(loc, 
             labels = c("index"), #can add multiple SNPs
             label_x = c(4, -5), 
             maxrows = 2, 
             filter_gene_biotype = 'protein_coding',
             highlight = Gene)
  dev.off()

}
