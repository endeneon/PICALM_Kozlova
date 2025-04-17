# Siwei 02 Jun 2024
# plot 5 western blot-confirmed genes at transcript level
# and their epitope sites;

# init #####
{
  library(Gviz)

  library(rtracklayer)
  library(GenomicFeatures)
  library(BSgenome)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(ensembldb)
  # library(org.Hs.eg.db)

  library(grDevices)
  library(gridExtra)

  library(RColorBrewer)

  library(readr)
  library(stringr)
}
#
options(ucscChromosomeNames = F)
options(Gviz.ucscUrl = "https://genome.ucsc.edu/cgi-bin/")

# make a txdb from gencode v35 #####
txdb_gencode_v35 <-
  makeTxDbFromGFF(file = "~/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf",
                  dataSource = "GENCODE",
                  organism = "hg38",
                  taxonomyId = 9606)
# organism(txdb_gencode_v35) <- "hg38"
# track_gencode_v35 <-
#   GeneRegionTrack()

track_gencode_v35 <-
  GeneRegionTrack("~/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf")
# ensdb_hg38 <-
#   EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

# ensdb_hs108 <-
#   makeTxDbFromEnsembl(organism = "Homo sapiens",
#                       release = 108)
columns(txdb_gencode_v35)
subset_txdb <-
  filter(txdb_gencode_v35,
         filter = ~ TXID == "ENST00000409469.*")

plot_data_grTrack <-
  exonsBy(txdb_gencode_v35,
          by = "gene",
          use.names = F)
# plot_data_grTrack <-
#   exonsBy(txdb_gencode_v35, "")
plot_data_grTrack <-
  unlist(plot_data_grTrack)
names(plot_data_grTrack)
elementMetadata(plot_data_grTrack)$transcript <-
  names(plot_data_grTrack)

# output_track <-
#   plot_data_grTrack[str_detect()]

## function ####
plot_gene_transcript <-
  function(chr, start, end, gene_id = "",
           mcols = 100, strand = "+",
           x_offset_1 = 0,
           x_offset_2 = 0,
           ylimit = 400,
           title_name = "") {
    cell_type <-
      GRanges(seqnames = Rle(chr),
              seqinfo = Seqinfo(seqnames = chr,
                                genome = "hg38"),
              ranges = IRanges(start = start,
                               end = end,
                               names = chr,
                               seqnames = chr),
              mcols = as.data.frame(mcols),
              strand = Rle(strand(strand)))


    iTrack <-
      IdeogramTrack(genome = genome(cell_type),
                    chromosome = as.character(unique(seqnames(cell_type))),
                    fontcolor = "black",
                    fontsize = 18)
    gTrack <-
      GenomeAxisTrack(col = "black",
                      fontcolor = "black",
                      fontsize = 16,
                      scale = 0.1)

    # mapped_gene_id <-
    #   mapIds(x = org.Hs.eg.db,
    #          keys = gene_name,
    #          column = )

    plot_df <-
      select(txdb_gencode_v35,
             keys = gene_id,
             columns = c("EXONCHROM",
                         "EXONSTART",
                         "EXONEND",
                         "EXONSTRAND",
                         "TXTYPE",
                         "GENEID",
                         "EXONID",
                         "TXNAME"),
             keytype = "GENEID")
    direct_plot_df <-
      plot_df
    colnames(direct_plot_df) <-
      c("gene",
        "exon",
        "chromosome",
        "strand",
        "start",
        "end",
        "transcript",
        "feature")
    grTrack <-
      GeneRegionTrack(direct_plot_df,
                      genome = "hg38",
                      chromosome = "chr1",
                      transcript = direct_plot_df$transcript)
    # z <-
    #   ranges(cell_type)
    # grTrack <-
    #   txdb_gencode_v35
    #   # track_gencode_v35
    # ranges(grTrack) <- z

    plotTracks(list(
      iTrack, gTrack, grTrack
      # gwasTrack,
      # # alTrack_hMG, alTrack_hAst,
      # alTrack_iMG, alTrack_iAst,
      # alTrack_NGN2,
      # alTrack_DN, alTrack_GA,
      # snpTrack,
      # grTrack
    ),
    # sizes = c(
    #   0.5, 0.5, 3
    #   # 1,
    #   # 1, 1,
    #   # # 1, 1,
    #   # 1,
    #   # 1,1,
    #   # 0.5,
    #   # 3
    # ),
    chromosome = cell_type@ranges@NAMES,
    from = (cell_type@ranges@start - x_offset_1),
    to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
    # transcriptAnnotation = "transcript",
    collapseTranscripts = "transcript",
    collapse = T,
    showID = T,
    geneSymbol = T,
    transcriptAnnotation = "transcript")#,

  }

plot_gene_transcript(chr = "chr7",
                     start = 148690000,
                     end = 148802000,
                     gene_id = "ENSG00000055130.17")

# add highlight #####
plot_gene_transcript <-
  function(chr, start, end,
           gene_name = "",
           gene_id = "",
           highlight_start = 0,
           highlight_width = 30,
           mcols = 100, strand = "+",
           x_offset_1 = 0,
           x_offset_2 = 0,
           ylimit = 400,
           edited_site = NA,
           title_name = "") {
    cell_type <-
      GRanges(seqnames = Rle(chr),
              seqinfo = Seqinfo(seqnames = chr,
                                genome = "hg38"),
              ranges = IRanges(start = start,
                               end = end,
                               names = chr,
                               seqnames = chr),
              mcols = as.data.frame(mcols),
              strand = Rle(strand(strand)))


    iTrack <-
      IdeogramTrack(genome = genome(cell_type),
                    chromosome = as.character(unique(seqnames(cell_type))),
                    fontcolor = "black",
                    fontsize = 18)
    gTrack <-
      GenomeAxisTrack(col = "black",
                      fontcolor = "black",
                      fontsize = 14)

    # mapped_gene_id <-
    #   mapIds(x = org.Hs.eg.db,
    #          keys = gene_name,
    #          column = )

    plot_df <-
      select(txdb_gencode_v35,
             keys = gene_id,
             columns = c("EXONCHROM",
                         "EXONSTART",
                         "EXONEND",
                         "EXONSTRAND",
                         "TXTYPE",
                         "GENEID",
                         "EXONID",
                         "TXNAME"),
             keytype = "GENEID")
    direct_plot_df <-
      plot_df
    colnames(direct_plot_df) <-
      c("gene",
        "exon",
        "chromosome",
        "strand",
        "start",
        "end",
        "transcript",
        "feature")
    print(direct_plot_df)
    grTrack <-
      GeneRegionTrack(direct_plot_df,
                      genome = "hg38",
                      chromosome = cell_type@ranges@NAMES,
                      name = "Transcripts",
                      # chromosome = "chr1",
                      transcript = direct_plot_df$transcript)

    htTrack <-
      HighlightTrack(trackList = list(iTrack,
                                      gTrack,
                                      grTrack),
                     chromosome = cell_type@ranges@NAMES,
                     start = c(highlight_start, edited_site),
                     width = c(highlight_width, 20),
                     genome = "hg38")
    # z <-
    #   ranges(cell_type)
    # grTrack <-
    #   txdb_gencode_v35
    #   # track_gencode_v35
    # ranges(grTrack) <- z

    plotTracks(list(
      htTrack),

    # chromosome = cell_type@ranges@NAMES,
    from = (cell_type@ranges@start - x_offset_1),
    to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
    # transcriptAnnotation = "transcript",
    collapseTranscripts = "transcript",
    collapse = T,
    showID = T,
    geneSymbol = T,
    transcriptAnnotation = "transcript",
    main = str_c(gene_name,
                 gene_id,
                 sep = ', '), cex.main = 1)#,

  }

plot_gene_transcript(chr = "chr8",
                     start = 52605000,
                     end = 52715000,
                     highlight_start = 52622400,
                     highlight_width = 100,
                     gene_name = "RB1CC1",
                     edited_site = 52683964,
                     gene_id = "ENSG00000023287.13")


plot_gene_transcript(chr = "chr6",
                     start = 156694000,
                     end = 157208000,
                     highlight_start = 157184400,
                     highlight_width = 50,
                     gene_name = "ARID1B",
                     edited_site = 156778374,
                     gene_id = "ENSG00000049618.24")

# plot_gene_transcript(chr = "chr6",
#                      start = 156776000,
#                      end = 156778422,
#                      highlight_start = 156776360,
#                      highlight_width = 50,
#                      gene_name = "ARID1B",
#                      gene_id = "ENSG00000271551.2")

plot_gene_transcript(chr = "chr14",
                     start = 21375000,
                     end = 21432000,
                     highlight_start = 21385416,
                     highlight_width = 50,
                     gene_name = "CHD8",
                     edited_site = 21429305,
                     gene_id = "ENSG00000100888.15")

plot_gene_transcript(chr = "chr7",
                     start = 21410000,
                     end = 21514900,
                     highlight_start = 21429000,
                     highlight_width = 1843,
                     gene_name = "SP4",
                     edited_site = 21429421,
                     gene_id = "ENSG00000105866.15")


plot_gene_transcript(chr = "chr7",
                     start = 148677000,
                     end = 148802000,
                     highlight_start = 148753800,
                     highlight_width = 12780,
                     gene_name = "CUL1",
                     edited_site = 148730185,
                     gene_id = "ENSG00000055130.17")

plot_gene_transcript(chr = "chr16",
                     start = 89235000,
                     end = 89500000,
                     highlight_start = 89281208,
                     highlight_width = 700,
                     edited_site = 89316947,
                     gene_name = "ANKRD11",
                     gene_id = "ENSG00000167522.16")

select(txdb_gencode_v35,
       keys = "ENSG00000055130.*",
       columns = c("EXONCHROM",
                   "EXONSTART",
                   "EXONEND",
                   "EXONSTRAND",
                   "TXTYPE",
                   "GENEID",
                   "EXONID",
                   "TXNAME"),
       keytype = "GENEID")

seqinfo(track_gencode_v35)

colnames(plot_data_grTrack)

columns(txdb_gencode_v35)
keys(txdb_gencode_v35,
     column = "")
genome(txdb_gencode_v35) <- 9606
plot_df <-
  select(txdb_gencode_v35,
         keys = "ENSG00000000938.13",
         columns = c("EXONCHROM",
                     "EXONSTART",
                     "EXONEND",
                     "EXONSTRAND",
                     "TXTYPE",
                     "GENEID",
                     "EXONID",
                     "TXNAME"),
         keytype = "GENEID")
colnames(plot_df)
range(txdb_gencode_v35)
gr_region_trk <-
  GeneRegionTrack(
    # start = min(plot_df$EXONSTART) - 5000,
    #               end = max(plot_df$EXONEND),
                  name = "test_gr",
                  rstarts = plot_df$EXONSTART,
                  rends = plot_df$EXONEND,
                  chromosome = "chr1",
                  genome = "hg38",
                  transcript = plot_df$TXID,
                  gene = plot_df$TXID,
                  strand = plot_df$EXONSTRAND,
                  symbol = plot_df$TXID,
                  exon = plot_df$EXONID,
                  feature = plot_df$TXID)
# identifier(gr_region_trk)
plotTracks(list(gr_region_trk),
           showID = T,
           collapse = F,
           # collapseTranscripts = "exon",
           # geneSymbol = T,
           transcriptAnnotation = "transcript")
direct_plot_df <-
  plot_df
colnames(direct_plot_df) <-
  c("gene",
    "exon",
    "chromosome",
    "strand",
    "start",
    "end",
    "transcript",
    "feature")
plotTracks(GeneRegionTrack(direct_plot_df,
                           genome = "hg38",
                           chromosome = "chr1",
                           transcript = direct_plot_df$transcript),
           showID = T,
           geneSymbol = T,
           transcriptAnnotation = "transcript")


data("cyp2b10")

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
keys(x = "CUL1",
          keytype = "SYMBOL")
AnnotationDbi::select(org.Hs.eg.db, keys = )
library(Hom)

