# Siwei 24 Apr 2024
# make Upset plots to show peaks between all cell types

# init #####
{
  library(readr)
  library(UpSetR)

  library(stringr)

  library(RColorBrewer)
}

# load data ####

df_peaks_file <-
  list.files(path = "indiv_peak_names",
             pattern = "*.bed",
             full.names = T)

df_raw <-
  vector(mode = "list",
         length = length(df_peaks_file))

for (i in 1:length(df_raw)) {
  print(df_peaks_file[i])
  df_peak_list_readin <-
    read.delim(file = df_peaks_file[i],
               header = F)
  df_raw[[i]] <-
    unlist(df_peak_list_readin)

}
names(df_raw) <-
  c("iAst",
    "iDN",
    "iGA",
    "iMG",
    "iGlut")

lst_data <-
  vector(mode = "list",
         length = length(df_raw))

# as.vector(unlist(df_raw[[1]]$`new iPSC line ID`))
for (i in 1:length(lst_data)) {
  lst_data[[i]] <- df_raw[[i]]
}

names(lst_data) <-
  c("iAst",
    "iDN",
    "iGA",
    "iMG",
    "iGlut")

upset(fromList(lst_data),
      sets = rev(c("iMG", "iAst", "iGA", "iDN", "iGlut")),
      keep.order = T,
      order.by = "degree",
      decreasing = F,
      sets.bar.color = "darkblue",
      sets.x.label = "Total peaks per cell type",
      mainbar.y.label = "Peaks in the category",
      set_size.angles = 45,
      number.angles = 0,
      shade.color = "black",
      queries = list(list(query = intersects,
                          params = list("iMG", "iAst"),
                          color = "darkred",
                          active = T),
                     # list(query = intersects,
                     #      params = list("iMG", "iAst", "iGlut"),
                     #      color = "darkred",
                     #      active = T),
                     list(query = intersects,
                          params = list("iGA", "iGlut", "iDN"),
                          color = "darkred",
                          active = T),
                     # list(query = intersects,
                     #      params = list("iDN"),
                     #      color = "black",
                     #      active = T),
                     list(query = intersects,
                          params = list("iMG", "iAst", "iGA", "iGlut", "iDN"),
                          color = "darkred",
                          active = T)
                     )
      )

upset(fromList(lst_data),
      sets = rev(c("iMG", "iAst", "iGA", "iDN", "iGlut")),
      keep.order = T,
      order.by = "freq",
      sets.bar.color = "darkblue",
      sets.x.label = "Cell line count",
      mainbar.y.label = "Lines in the category",
      queries = list(list(query = intersects,
                          params = list("iMG", "iAst", "iGA"),
                          color = "darkred",
                          active = T)
      )
)
###
test_list <-
  list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13),
       two = c(1, 2, 4, 5, 10),
       three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

upset(movies,
      set.metadata = list(data = metadata,
                          plots = list(list(type = "hist",
                                            column = "avgRottenTomatoesScore",
                                            assign = 20),
                                       list(type = "bool",
                                            column = "accepted",
                                            assign = 5,
                                            colors = c("#FF3333", "#006400")),
                                       list(type = "text", column = "Cities",
                                            assign = 5,
                                            colors = c(Boston = "green",
                                                       NYC = "navy",
                                                       LA = "purple")),
                                       list(type = "matrix_rows",
                                            column = "Cities",
                                            colors = c(Boston = "green",
                                                       NYC = "navy", LA = "purple"),
                                            alpha = 0.5))),
      queries = list(list(query = intersects,
                          params = list("Drama"),
                          color = "red",
                          active = F),
                     list(query = intersects,
                          params = list("Action", "Drama"),
                          active = T),
                     list(query = intersects,
                          params = list("Drama", "Comedy", "Action"),
                          color = "orange",
                          active = T)),
      attribute.plots = list(gridrows = 45,
                             plots = list(list(plot = scatter_plot,
                                               x = "ReleaseDate",
                                               y = "AvgRating",
                                               queries = T),
                                          list(plot = scatter_plot,
                                               x = "AvgRating",
                                               y = "Watches",
                                               queries = F)),
                             ncols = 2),
      query.legend = "bottom")
