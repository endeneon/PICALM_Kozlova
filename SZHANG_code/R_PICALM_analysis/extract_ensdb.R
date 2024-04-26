library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
columns(ensdb)
head(keys(ensdb))
head(keytypes(ensdb))
keytypes(ensdb)
head(keys(ensdb,keytype = "GENEBIOTYPE"))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <-
  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
keytypes(txdb)
columns(txdb)

head(txdb)
col_info <-
  columns(txdb)
keys(txdb)
head(transcripts(txdb))

AnnotationDbi::select(txdb,
                      keys = as.character(unlist(head(transcripts(txdb)$tx_id))),
                      columns = c("TXID", "TXTYPE", "GENEID"),
                      keytype = "TXID")

library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
columns(ensdb)
keytypes(ensdb)

AnnotationDbi::select(ensdb,
                      keys = as.character(unlist(head(transcripts(ensdb)$tx_name))), # entry
                      columns = c("TXID", "TXBIOTYPE", "GENEID", "TXNAME"), # return values
                      keytype = "TXID") # lookup word

full_data_df <-
  AnnotationDbi::select(ensdb,
                        keys = as.character(unlist(transcripts(ensdb)$tx_name)), # entry
                        columns = c("TXID", "TXBIOTYPE", "GENEID"), # return values
                        keytype = "TXID") # lookup word
unique(full_data_df$TXBIOTYPE)

protein_tx_df <-
  full_data_df[full_data_df$TXBIOTYPE == "protein_coding", ]

head(transcripts(ensdb))

library(tximport)
tx_lookup_table <-
  protein_tx_df[, c(1, 3)]
colnames(tx_lookup_table)
tximport(tx2gene = tx_lookup_table,
         txIdCol = "TXID",
         geneIdCol = "GeneID", 
         ignoreTxVersion = T)
