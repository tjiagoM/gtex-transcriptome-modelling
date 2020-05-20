#Load the library
library(sva)

# Run script as:
# Rscript --no-save --no-restore --verbose 02_01_correct_confounds.R > output_02_01.txt 2>&1

#GENES IN THE ROWS AND SAMPLES IN THE COLUMNS
tissues <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala', 
             'Brain_Anterior_cingulate', 'Brain_Caudate', 'Brain_Cerebellar', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex', 'Brain_Hippocampus', 
             'Brain_Hypothalamus', 'Brain_Nucleus', 'Brain_Putamen', 'Brain_Spinal_cord', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Cultured', 
             'Cells_EBV', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastro', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial', 'Heart_L_Vent', 
             'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Epsd', 
             'Skin_Sun_Epsd', 'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood')

for (tissue in tissues)
{
  print(paste('Tissue ', tissue, ': ', sep=''))
  # row.names = 1 specifies that 1st column is the column with the rownames
  file=read.csv(file=paste('data_filtered/only_geneids_', tissue, '.csv', sep=''), row.names = 1)
  # Create the model
  mod=matrix(1,nrow=dim(file)[1],ncol=1)
  colnames(mod)="Intercept"
  #Set the rowname of the file as the first column of the matrix of gene expressions 
  n.pc=num.sv(t(file),mod, method = "be")
  print(paste(n.pc, 'components'))
  #Create the new adjusted dataset
  dat.adjusted = sva_network(file, n.pc)
  #Write the new dataframe with the metrics adjusted into a csv file 
  write.csv(dat.adjusted,file=paste('data_filtered/only_geneids_CORRECTED_', tissue, '.csv', sep=''))
} 
