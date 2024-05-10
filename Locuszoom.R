library(locuszoomr); library(EnsDb.Hsapiens.v75); library(data.table)

#Change input for your analysis

Pop="EUR"
Gene="PAX8"
Assay="IL1RN_P18510_OID20700_v1_Inflammation"

#make object for Filename.Input and Filename.Output
Filename.Input=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/input/results_", Gene, "_region_", Assay, "_rsID.txt")
Filename.Output=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/output/locuszoom_", Gene, "_region_", Assay, "_", Pop, ".pdf")

#########################################################################

#read input 
dat <- fread(Filename.Input)

#create locus object for plotting
#gene region +/-500kb 
loc <- locus(gene = Gene, data=dat, flank = 5e5, ens_db="EnsDb.Hsapiens.v75")

#add LD: to do this, get the token in here: https://ldlink.nih.gov/?tab=apiaccess
#it queries LDlink (https://ldlink.nci.nih.gov/) via the LDlinkR package.
#pop: a 1000 Genomes Project population. multiple populations allowed.
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

###########################################
#For loop 
Assay=c("IL1RN_P18510_OID20700_v1_Inflammation", "IL1R1_P14778_OID21116_v1_Neurology", "HEPACAM2_A8MVW5_OID31151_v1_Neurology_II", "CALCA_P01258_OID20983_v1_Neurology", "IL36G_Q9NZH8_OID30551_v1_Inflammation_II", "MMP1_P03956_OID20672_v1_Inflammation", "IFIT3_O14879_OID31283_v1_Oncology_II", "KLK15_Q9H2R5_OID31427_v1_Oncology_II", "TG_P01266_OID31359_v1_Oncology_II", "DPP10_Q8N608_OID20527_v1_Inflammation", "UMOD_P07911_OID20237_v1_Cardiometabolic")

for (i in 1:length(Assay)) {

  #make object for Filename.Input and Filename.Output
  Filename.Input=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/input/results_", Gene, "_region_", Assay[i], "_rsID.txt")
  Filename.Output=paste0("/mnt/zfs/P201_UKB_pQTL/results/", Gene, "/region/plot/output/locuszoom_", Gene, "_region_", Assay[i], "_", Pop, ".pdf")
  
  #read input 
  dat <- fread(Filename.Input)
  
  #create locus object for plotting
  #gene region +/-500kb 
  loc <- locus(gene = Gene, data=dat, flank = 5e5, ens_db="EnsDb.Hsapiens.v75")
  
  #add LD: to do this, get the token in here: https://ldlink.nih.gov/?tab=apiaccess
  #it queries LDlink (https://ldlink.nci.nih.gov/) via the LDlinkR package.
  #pop: a 1000 Genomes Project population. multiple populations allowed.
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
