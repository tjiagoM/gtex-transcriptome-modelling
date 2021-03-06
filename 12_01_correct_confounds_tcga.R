#Load the library
library(sva)

# Run script as:
# Rscript --no-save --no-restore --verbose 12_01_correct_confounds_tcga.R > outputs/output_12_01.txt 2>&1

tissues <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')


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
