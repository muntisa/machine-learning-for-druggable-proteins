library(Rcpi)

# read S2SNet files
df <- read.table(text = readLines("seqs.Screening_2_OncoOmics_Genes.txt", warn = FALSE), sep = "\t")

# get sequences as vector
seqs<-apply(df['V3'], 1, as.vector)
# no of seqs BEFORE protein check
ini <- length(seqs)

# To assure that the protein sequences only have the twenty standard amino acid types
# which is required for the descriptor computation, we use the checkProt() function
# to do the amino acid type sanity checking and remove the non-standard protein sequences:
seqs <- seqs[(sapply(seqs, checkProt))]
# no of seqs AFTEER protein check
fin <- length(seqs)

# deleted sequences
fin - ini

### Tripeptide composition Descriptor 
#################################################### 
# extractProtTC
desc_TC = t(sapply(seqs, extractProtTC))
# scale
#TC_s<-scale(desc_TC)
# convert to dataframe
TC_s<-as.data.frame(desc_TC)
# copy rownames as column
TC_s <- cbind(rownames(TC_s), data.frame(TC_s, row.names=NULL))
# change name of 1st column
colnames(TC_s)[1] = "V3"
# merge 2 df using a col
TC_lab<-merge(TC_s, df)
# write the result as file
write.csv(TC_lab, file = "TC_seqs.Screening_2_OncoOmics_Genes.csv")
