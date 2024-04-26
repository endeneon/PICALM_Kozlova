# Siwei 18 Feb 2024
# make a GRangesObject of TSS -2/+1 kb
# Need to test which genes to include (if only protein_coding, 19 k)

# init ####
{
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(org.Hs.eg.db)



  library(stringr)
  library(future)

  # library(MASS)

}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# funcs ####
CollapseToLongestTranscript_2 <-
  function(ranges) {
    range.df <- data.table::as.data.table(x = ranges)
    range.df$strand <- as.character(x = range.df$strand)
    range.df$strand <- ifelse(test = range.df$strand == "*",
                              yes = "+", no = range.df$strand)
    collapsed <- range.df[, .(unique(seqnames), min(start),
                              max(end), strand[[1]], gene_biotype[[1]], gene_name[[1]]),
                          "gene_id"]
    colnames(x = collapsed) <- c("gene_id", "seqnames", "start",
                                 "end", "strand", "gene_biotype", "gene_name")
    collapsed$gene_name <- .Internal(make.unique(names = collapsed$gene_name,
                                                 sep = "."))
    gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = collapsed,
                                                           keep.extra.columns = TRUE)
    return(gene.ranges)
  }


# load a hg38 annot genome
hg38_annot <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(hg38_annot) <-
  paste0('chr',
         seqlevels(hg38_annot))
genome(hg38_annot) <-
  'hg38'

# keep autosomes and sort by chr names
hg38_annot <-
  keepStandardChromosomes(hg38_annot,
                          pruning.mode = "coarse")
hg38_annot <-
  sortSeqlevels(hg38_annot)


# check how many gene types in hg38_annot
unique(hg38_annot@elementMetadata@listData$gene_biotype)
unique(hg38_annot@elementMetadata@listData$type)
length(hg38_annot@elementMetadata@listData$gene_biotype)
length(unique(hg38_annot@elementMetadata@listData$gene_id))

# remove rRNAs
hg38_annot <-
  hg38_annot[hg38_annot@elementMetadata@listData$gene_biotype %in% c("protein_coding",
                                                                     "lincRNA",
                                                                     "processed_transcript"), ]
length(hg38_annot)

ens_hg38_transcript <-
  CollapseToLongestTranscript_2(hg38_annot)


columns(EnsDb.Hsapiens.v86)
head(select(EnsDb.Hsapiens.v86,
            keys = head(keys(EnsDb.Hsapiens.v86,
                             keytype = "GENEID")),
            columns = c("TXNAME",
                        # "TXTYPE",
                        "GENENAME"),
            keytype = "GENEID"))
# hg38_ensdb_annot_promoter <-


## use TxDb #####
hg38_txdb <-
  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
hg38_txdb <-
  EnsDb.Hsapiens.v86
columns(hg38_txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"   "CDSSTRAND"  "EXONCHROM"  "EXONEND"
# [10] "EXONID"     "EXONNAME"   "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"
# [19] "TXNAME"     "TXSTART"    "TXSTRAND"   "TXTYPE"
unique(hg38_txdb)
keys(hg38_txdb)
keys(hg38_txdb, keytype = "TXTYPE")
keytypes(hg38_txdb)
# [1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"     "TXNAME"

head(keys(hg38_txdb,
          keytype = "TXNAME"))

# keep autosomes and sort by chr names
hg38_txdb <-
  keepStandardChromosomes(hg38_txdb,
                          pruning.mode = "coarse")
hg38_txdb <-
  sortSeqlevels(hg38_txdb)

length(promoters(hg38_txdb,
               upstream = 2000,
               downstream = 1000,
               use.names = T))
head(select(hg38_txdb,
            keys = head(keys(hg38_txdb,
                             keytype = "GENEID")),
            columns = c("TXNAME",
                        "TXTYPE",
                        "GENEID"),
            keytype = "GENEID"))
unique(select(hg38_txdb,
            keys = keys(hg38_txdb,
                             keytype = "GENEID"),
            columns = c("TXNAME",
                        "GENEID"),
            keytype = "GENEID")$GENEID)
# length(select(hg))


## use org.Hs.eg.db
hg38_orghsegdb <-
  org.Hs.eg.db
columns(hg38_orghsegdb)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"
# [9] "EVIDENCEALL"  "GENENAME"     "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"
# [17] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"
# [25] "UCSCKG"       "UNIPROT"
keys(hg38_orghsegdb)
keytypes(hg38_orghsegdb)

unique(select(hg38_orghsegdb,
              keys = keys(hg38_orghsegdb,
                          keytype = "GENETYPE"),
              columns = c("ENSEMBL",
                          "GENETYPE"),
              keytype = "GENETYPE")$GENETYPE)
head(select(hg38_orghsegdb,
              keys = keys(hg38_orghsegdb,
                          keytype = "ENSEMBLTRANS"),
              columns = c("ENSEMBLTRANS",
                          "GENETYPE"),
              keytype = "ENSEMBLTRANS"))

protein_ncRNA <-
  select(hg38_orghsegdb,
              keys = keys(hg38_orghsegdb,
                          keytype = "ENSEMBLTRANS"),
              columns = c("ENSEMBLTRANS",
                          "GENETYPE"),
              keytype = "ENSEMBLTRANS")
protein_ncRNA <-
  protein_ncRNA[protein_ncRNA$GENETYPE %in% c("protein-coding",
                                                "ncRNA"), ]


## extract hg38 TX promoters and select
# ens_hg38_transcript <-
#   CollapseToLongestTranscript_2(hg38_annot)

# hg38_tx_promoters <-
#   promoters(hg38_txdb,
#             upstream = 2000,
#             downstream = 1000,
#             use.names = T)
hg38_tx_promoters <-
  promoters(ens_hg38_transcript,
            upstream = 2000,
            downstream = 1000,
            use.names = T)
unique(hg38_tx_promoters@elementMetadata$gene_biotype)
# remove the "."
# hg38_tx_promoters@ranges@NAMES <-
#   str_split(string = hg38_tx_promoters@ranges@NAMES,
#             pattern = "\\.",
#             simplify = T)[, 1]
# hg38_tx_promoters@ranges@NAMES

# hg38_tx_promoters@elementMetadata@listData$tx_name <-
#   str_split(string = hg38_tx_promoters@elementMetadata@listData$tx_name,
#             pattern = "\\.",
#             simplify = T)[, 1]
hg38_tx_promoters@elementMetadata@listData$tx_name <-
  hg38_tx_promoters@elementMetadata@listData$gene_name
hg38_tx_promoters@elementMetadata@listData$tx_id <-
  hg38_tx_promoters@elementMetadata@listData$tx_name
length(unique(hg38_tx_promoters@elementMetadata@listData$tx_name))

hg38_tx_reduced_promoters <-
  GenomicRanges::reduce(hg38_tx_promoters,
                        drop.empty.ranges = T,
                        # na.rm = T,
                        with.revmap = T,
                        ignore.strand = T,
                        min.gapwidth = 100)
revamp_index <-
  hg38_tx_reduced_promoters$revmap@partitioning@end

hg38_tx_reduced_promoters$revmap_annot  <-
  hg38_tx_promoters$tx_name[hg38_tx_reduced_promoters$revmap@partitioning@end]
save(hg38_tx_reduced_promoters,
     file = "hg38_tx_reduced_TSS_region_2k_1k.RData")
# filter
# hg38_tx_promoters <-
#   hg38_tx_promoters[hg38_tx_promoters@ranges@NAMES %in% protein_ncRNA$ENSEMBLTRANS, ]
# length(hg38_tx_promoters)
# coverage(hg38_tx_promoters)
# hg38_tx_promoters@ranges@width
# GenomicRanges::coverage(hg38_tx_promoters)
