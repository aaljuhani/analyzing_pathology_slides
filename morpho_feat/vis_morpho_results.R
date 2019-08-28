require(dplyr)    
library(plyr)
library(ggpubr)


setwd("~/Documents/projects/histomicstk_features/analyzing_pathology_slides/morpho_feat")

file_name = "Case 2_001.svs_40_morpho.csv"
morpho <- read.table(file_name,
                       header=TRUE, sep=",")

morpho$id <- seq.int(nrow(morpho))

morpho_df <- tbl_df(morpho)



morpho_df$tall_cell <- ifelse(morpho_df$Size.MajorAxisLength > 3* morpho_df$Size.MinorAxisLength , 3,
                              ifelse(morpho_df$Size.MajorAxisLength > 2* morpho_df$Size.MinorAxisLength, 2, 1))

morpho_df$tall_cell = factor(morpho_df$tall_cell)
summary(morpho_df)
str(morpho_df)
table(morpho_df$tall_cell)
table(morpho_df$tile, morpho_df$tall_cell)
tallcell_per_tile <-as.data.frame(table(morpho_df$tile, morpho_df$tall_cell))
str(tallcell_per_tile)


sort(morpho_df$Size.Area , decreasing = TRUE)


#summary(twice_tallcells)

#par(mfrow=c(2,3))


ggdensity(morpho_df, x = "Shape.Eccentricity",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          ylim = c(0, 30),
          title = file_name)




ggdensity(morpho_df, x = "Size.Perimeter",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggdensity(morpho_df, x = "Size.Area",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggdensity(morpho_df, x = "Shape.Circularity",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggdensity(morpho_df, x = "Shape.Eccentricity",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          ylim = c(0, 30),
          title = "Case 5")

ggdensity(morpho_df, x = "Shape.EquivalentDiameter",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggdensity(morpho_df, x = "Shape.Solidity",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggdensity(morpho_df, x = "Size.MajorAxisLength",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))


ggdensity(morpho_df, x = "Shape.MinorMajorAxisRatio",
          add = "mean", rug = TRUE,
          color = "tall_cell", fill = "tall_cell",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))


ggdensity(morpho_df, x = "tall_cell",
          add = "mean", rug = TRUE,
          color = "tile", fill = "tile")


ggplot(tallcell_per_tile, aes(x = "Var1", y = "Var1", fill = "Freq")) + 
  geom_bar(stat = "identity") +
  xlab("\nTiles") +
  ylab("Tall Cells\n") +
  guides(fill = FALSE) +
  theme_bw()

plot(tallcell_per_tile)

ggplot(data= tallcell_per_tile, aes(x = Var1, y= Freq, fill= Var2))+ 
  geom_bar(stat = "identity")+
  xlab("\nTiles") +
  ylab("Tall Cells\n") +
  theme_bw()

gghistogram(tallcell_per_tile, x = "Var1", y = "Freq", 
         add = "mean", color = "tall_cell", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#FF0000"),
         ylab = "Area", xlab = "Tiles")







hist(morpho_df$tall_cell , main= "Tall Cells in all tiles", xlab="area", col=morpho_df$tall_cell)



#plot all cells - Area
plot( morpho_df$tall_cell , morpho_df$Size.Area , main= "Cells in all tiles", xlab="area", col=c('gray', 'aquamarine2', 'red'))

# plot all cells - Circularity
barplot(sort(morpho_df$Shape.Circularity), main= "Cells in all tiles", xlab="Shape.Circularity")

# MajorAxisLength
barplot(sort(morpho_df$Size.MajorAxisLength), main= "Cells in all tiles", xlab="MajorAxisLength")

# MinorAxisLength
barplot(sort(morpho_df$Size.MinorAxisLength), main= "Cells in all tiles", xlab="MinorAxisLength")

# Perimeter
barplot(sort(morpho_df$Size.Perimeter), main= "Cells in all tiles", xlab="Perimeter")



# plot number of cells per tile 
groupby_tile <- group_by(morpho_df, tile)
cells_per_tile <- summarise(groupby_tile, cells = n_distinct(X))
cells_per_tile[order(cells_per_tile$cells)]
barplot(sort(cells_per_tile$cells), main="Number of Cells in Tiles", xlab="tile")



