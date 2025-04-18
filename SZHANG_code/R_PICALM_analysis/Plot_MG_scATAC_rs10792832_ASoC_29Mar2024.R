# Siwei 29 Mar 2024
# Plot CLU1 SNP rs1532278

# init #####
{
  library(Gviz)

  library(rtracklayer)
  # library(BSgenome)
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
options(Gviz.ucscUrl = "https://genome.ucsc.edu/cgi-bin/hgGateway")


plot_AsoC_peaks <-
  function(chr, start, end, gene_name = "",
                            alTrackRegion,
                            alTrackRefT,
                            alTrackAltC,
                            mcols = 100, strand = "+",
                            x_offset_1 = 0, x_offset_2 = 0, ylimit = 400,
                            title_name = "") {
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

  # iTrack <- IdeogramTrack(genome = genome(cell_type),
  #                         chromosome = as.character(unique(seqnames(cell_type))),
  #                         fontcolor = "black",
  #                         fontsize = 18)

  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16,
                            scale = 0.4)

  alTrack_Region <-
    AlignmentsTrack(alTrackRegion,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/subsampled_0hr_GABA.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "palegreen",
                    fill.coverage = "palegreen",
                    col.coverage = "palegreen",
                    ylim = c(0, ylimit),
                    cex.axis = 0.75,
                    col.title = "black",
                    col.axis = "black",
                    # transformation = function(x) {x / 15},
                    name = title_name,
                    fontsize = 12,
                    show.title = FALSE)

  alTrack_T <-
    AlignmentsTrack(alTrackRefT,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.ref.T.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "darkred",
                    fill.coverage = "darkred",
                    col.coverage = "darkred",
                    ylim = c(0, ylimit),
                    cex.axis = 0.75,
                    # transformation = function(x) {x / 15},
                    # col.title="black",
                    name = "Allele_T")

  alTrack_C <- AlignmentsTrack(alTrackAltC,
                               # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.alt.C.bam",
                               isPaired = F, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkblue",
                               fill.coverage = "darkblue",
                               col.coverage = "darkblue",
                               ylim = c(0, ylimit),
                               cex.axis = 0.75,
                               # transformation = function(x) {x / 15},
                               # col.title="black",
                               name = "Allele_C")

  OvlTrack <- OverlayTrack(trackList = list(alTrack_Region,
                                            alTrack_T,
                                            alTrack_C),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  snpTrack <- AnnotationTrack(start = 86156833,
                              end = 86156833, chromosome = "chr11",
                              id = "rs10792832", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs10792832"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              # fontsize = 14,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")


  plotTracks(list(gTrack,
                  OvlTrack,
                  snpTrack),
             sizes = c(0.5,
                       1.5,
                       0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}

# plot_AsoC_peaks(chr = "chr8",
#                 start = 9787096,
#                 end = 9789096,
#                 alTrackRegion = TrackRegionBam[i],
#                 alTrackRefT = TrackRefTBam[i],
#                 alTrackAltC = TrackAltCBam[i],
#                 ylimit = 100,
#                 title_name = BAM_names[i])

plot_AsoC_peaks(chr = "chr11",
                start = 86156733,
                end = 86156933,
                x_offset_1 = 30,
                x_offset_2 = 30,
                alTrackRegion = "/data/FASTQ/sage_synapse/all_het_rs10792832_merged_sorted.bam",
                alTrackRefT = "/data/FASTQ/sage_synapse/split_rs10792832/.ref.A.bam",
                alTrackAltC = "/data/FASTQ/sage_synapse/split_rs10792832/.alt.G.bam",
                ylimit = 200,
                title_name = "snATAC")
