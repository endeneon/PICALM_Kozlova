# Siwei 25 Jun 2024
# plot rs2027349 using the bw files from
# https://personal.broadinstitute.org/bjames/AD_snATAC/bigWig_TSS6/
# and their epitope sites;

# init #####
{
  library(Gviz)

  library(rtracklayer)
  library(GenomicFeatures)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ensembldb)
  library(org.Hs.eg.db)

  library(grDevices)
  library(gridExtra)

  library(RColorBrewer)

  library(readr)
  library(stringr)
}
#

options(Gviz.ucscUrl = "https://genome.ucsc.edu/cgi-bin/")
options(ucscChromosomeNames = F)

# make a txdb from gencode v35 #####
txdb_gencode_v35 <-
  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene


# make the master function ####

plot_AD_snATAC <-
  function(chr, start, end,
           gene_name = "",
           SNP_position = 0,
           SNP_ID = "",
           mcols = 100, strand = "+",
           x_offset_1 = 0, x_offset_2 = 0,
           track_type = "l",
           font_size = 10,
           show_axis = F,
           iMG_ylimit = 800,
           ylimit = 800) {
  cell_type <-
    GRanges(seqnames = Rle(chr),
            seqinfo = Seqinfo(seqnames = chr,
                              genome = "hg38"),
            ranges = IRanges(start = start,
                             end = end,
                             names = chr),
            mcols = as.data.frame(mcols),
            strand = Rle(strand(strand)))

  print(as.character(unique(seqnames(cell_type))))

  iTrack <- IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))),
                          fontcolor = "black",
                          fontsize = 16)
  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 10
                            # scale = 0.1 # not using scale, draw full axis
  )

  alTrack_iMG <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation = function(x) {x * 5},
                    ylim = c(0, iMG_ylimit),
                    fontcolor = "black",
                    fontcolor.title = "black",
                    background.title = brewer.pal(n = 9,
                                                  name = "Set1")[2],
                    fill.coverage = brewer.pal(n = 9,
                                               name = "Set1")[2],
                    col.coverage = brewer.pal(n = 9,
                                              name = "Set1")[2],
                    fontsize = font_size,
                    col.axis = "transparent",
                    showAxis = F,
                    name = "iMG")

  dtTrack_Ast <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/Ast_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 9,
                                            name = "Set1")[1],
              fill = brewer.pal(n = 9,
                                         name = "Set1")[1],
              col = brewer.pal(n = 9,
                                name = "Set1")[1],
              col.histogram = brewer.pal(n = 9,
                               name = "Set1")[1],
              fill.coverage = brewer.pal(n = 9,
                                         name = "Set1")[1],
              col.coverage = brewer.pal(n = 9,
                                        name = "Set1")[1],
              fill.horizon = brewer.pal(n = 9,
                                        name = "Set1")[1],
              fill.histogram = brewer.pal(n = 9,
                                          name = "Set1")[1],
              fontcolor = "black",
              fontcolor.title = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 2},
              # cex.axis = 2,
              fontsize = font_size,
              type = track_type,
              showAxis = show_axis,
              name = "Ast")

  dtTrack_Ex <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/Ex_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 9,
                                            name = "Set1")[3],
              fill.coverage = brewer.pal(n = 9,
                                         name = "Set1")[3],
              fill.horizon = brewer.pal(n = 9,
                                        name = "Set1")[3],
              fill.histogram = brewer.pal(n = 9,
                                        name = "Set1")[3],
              col.coverage = brewer.pal(n = 9,
                                        name = "Set1")[3],
              col.histogram = brewer.pal(n = 9,
                                        name = "Set1")[3],
              fontcolor = "black",
              ylim = c(0, ylimit),
              fontcolor.title = "black",
              # transformation = function(x) {x * 2},
              # cex.axis = 2,
              fontsize = font_size,
              type = track_type,
              showAxis = show_axis,
              name = "Ex")

  dtTrack_In <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/In_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 8,
                                            name = "Paired")[3],
              fill.coverage = brewer.pal(n = 8,
                                         name = "Paired")[3],
              fill.horizon = brewer.pal(n = 8,
                                         name = "Paired")[3],
              fill.histogram = brewer.pal(n = 8,
                                        name = "Paired")[3],
              col.coverage = brewer.pal(n = 8,
                                        name = "Paired")[3],
              col.histogram = brewer.pal(n = 8,
                                        name = "Paired")[3],
              fontcolor = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 2},
              fontcolor.title = "black",
              # cex.axis = 2,
              fontsize = font_size,
              type = track_type,
              showAxis = show_axis,
              name = "In")

  dtTrack_MG <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/Microglia_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 10,
                                            name = "Paired")[9],
              fill.coverage = brewer.pal(n = 10,
                                         name = "Paired")[9],
              fill.horizon = brewer.pal(n = 10,
                                         name = "Paired")[9],
              fill.histogram = brewer.pal(n = 10,
                                        name = "Paired")[9],
              col.coverage = brewer.pal(n = 10,
                                        name = "Paired")[9],
              col.histogram = brewer.pal(n = 10,
                                        name = "Paired")[9],
              fontcolor = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 8},
              fontcolor.title = "black",
              # cex.axis = 2,
              fontsize = font_size,
              showAxis = show_axis,
              type = track_type,
              name = "Microglia")

  dtTrack_Olig <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/Oligo_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 10,
                                            name = "Paired")[8],
              fill.coverage = brewer.pal(n = 10,
                                         name = "Paired")[8],
              fill.horizon = brewer.pal(n = 10,
                                         name = "Paired")[8],
              fill.histogram = brewer.pal(n = 10,
                                        name = "Paired")[8],
              col.coverage = brewer.pal(n = 10,
                                        name = "Paired")[8],
              col.histogram = brewer.pal(n = 10,
                                        name = "Paired")[8],
              fontcolor = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 2},
              fontcolor.title = "black",
              # cex.axis = 2,
              fontsize = font_size,
              showAxis = show_axis,
              type = track_type,
              name = "Oligo")

  dtTrack_OPC <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/OPC_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 9,
                                            name = "Set1")[7],
              fill.coverage = brewer.pal(n = 9,
                                         name = "Set1")[7],
              fill.horizon = brewer.pal(n = 9,
                                         name = "Set1")[7],
              fill.histogram = brewer.pal(n = 9,
                                        name = "Set1")[7],
              col.coverage = brewer.pal(n = 9,
                                        name = "Set1")[7],
              col.histogram = brewer.pal(n = 9,
                                        name = "Set1")[7],
              fontcolor = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 2},
              fontcolor.title = "black",
              # cex.axis = 2,
              fontsize = font_size,
              showAxis = show_axis,
              type = track_type,
              name = "OPC")

  dtTrack_Vas <-
    DataTrack(range = "AD_snATAC_bigWig_TSS6/PerEndo_nonAD-TileSize-50-normMethod-nFrags-ArchR.bw",
              chromosome = as.character(unique(seqnames(cell_type))),
              genome = "hg38",
              stream = T,
              background.title = brewer.pal(n = 9,
                                            name = "Set1")[9],
              fill.coverage = brewer.pal(n = 9,
                                         name = "Set1")[9],
              fill.horizon = brewer.pal(n = 9,
                                         name = "Set1")[9],
              fill.histogram = brewer.pal(n = 9,
                                        name = "Set1")[9],
              col.coverage = brewer.pal(n = 9,
                                        name = "Set1")[9],
              col.histogram = brewer.pal(n = 9,
                                        name = "Set1")[9],
              fontcolor = "black",
              ylim = c(0, ylimit),
              # transformation = function(x) {x * 2},
              fontcolor.title = "black",
              # cex.axis = 2,
              fontsize = font_size,
              showAxis = show_axis,
              type = track_type,
              name = "Vas")


  # snpTrack <- AnnotationTrack(start = c(86065502, 86156833), # hg19 = c("85867875)
  #                             end = c(86065502, 86156833),
  #                             chromosome = "chr11", # rs78710909
  #                             id = c("rs867611", "rs10792832"),
  #                             shape = "box",
  #                             name = "SNP", strand = "*",
  #                             group = c("rs867611", "rs10792832"),
  #                             fontcolor.group = "black", fontcolor.item = "black",
  #                             fontsize = 16,
  #                             fontsize.group = 10,
  #                             col = "black", col.title = "black",
  #                             just.group = "below",
  #                             showID = TRUE,
  #                             cex.group = 1,
  #                             groupAnnotation = "id")

  snpTrack <- AnnotationTrack(start = SNP_position, # hg19 = c("85867875)
                              end = SNP_position,
                              chromosome = chr, # rs78710909
                              id = SNP_ID,
                              shape = "box",
                              name = "SNP", strand = "*",
                              group = SNP_ID,
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              fontsize.group = 10,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")

  ########### plotting
  ucscGenes <- UcscTrack(genome=genome(cell_type), table="ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                             mcols(z)$symbol), "SYMBOL","REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor <- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 16

  ######

  # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
  #                           start = c(26212000, 26241000),
  #                           width = c(15000, 2000)
  #                           chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(iTrack, gTrack,
                  alTrack_iMG,
                  dtTrack_Ast, dtTrack_Ex, dtTrack_In, dtTrack_MG,
                  dtTrack_Olig, dtTrack_OPC, dtTrack_Vas,
                  # alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS,
                  snpTrack, grTrack),
             sizes = c(0.5,0.75,
                       1,
                       1,1,1,1,1,1,1,
                       0.5,1),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}

plot_AD_snATAC(chr = "chr11",
               start = 85950000,
               end = 86200000,
               SNP_position = 86156833,
               SNP_ID = "rs10792832",
               track_type = "hist",
               font_size = 16,
               iMG_ylimit = 350,
               ylimit = .05)

plot_AD_snATAC(chr = "chr8",
               start = 27603798,
               end = 27609398,
               SNP_position = 27608898,
               SNP_ID = "rs1532278",
               track_type = "hist",
               font_size = 16,
               iMG_ylimit = 350,
               ylimit = .05)

plot_AD_snATAC(chr = "chr8",
               start = 27583329,
               end = 27623318,
               SNP_position = 27608898,
               SNP_ID = "rs1532278",
               track_type = "hist",
               font_size = 16,
               iMG_ylimit = 350,
               ylimit = .05)

plot_AD_snATAC(chr = "chr8",
               start = 27583329,
               end = 27623318,
               SNP_position = 27608898,
               SNP_ID = "rs1532278",
               track_type = "hist",
               font_size = 16,
               iMG_ylimit = 350,
               ylimit = .05)

plot_AD_snATAC(chr = "chr5",
               start = 86600000,
               end = 87050000,
               SNP_position = 86997803,
               SNP_ID = "rs62375391",
               track_type = "hist",
               font_size = 16,
               iMG_ylimit = 350,
               ylimit = .05)
