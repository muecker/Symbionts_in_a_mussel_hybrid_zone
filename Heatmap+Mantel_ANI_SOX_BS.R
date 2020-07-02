# Plot heatmap of Bathymodiolus SOX symbiont ANI values and test correlation with sampling year & host species 2020-02-05
# Merle Ücker

# Load libraries
library(gplots)
library(RColorBrewer)
library(vegan)

# Load data
ani <- read.table("Average_nucleotide_identity_SOX_symbionts_Broken_Spur.csv", sep = "\t", header = TRUE, check.names = FALSE)

# Prepare matrix
rownames <- ani[,1] # set row names to MAG names
row.names(ani) <- rownames
ani[,1] <- NULL # remove names column
ani[is.na(ani)] <- 0 # set missing data to 0
ani_matrix <- as.matrix(ani)

# Load metadata
metadata <- read.table("metadata2", sep="\t", header = TRUE)
# metadata2:
# Species	Year
# Bhyb	2001
# Bhyb	2001
# Bput	2001
# Bhyb	1997
# Bput	2001
# Bhyb	1997
# Bput	1997
# Bput	2001
# Bhyb	1997
# Bhyb	1997
# Bhyb	1997
# Bput	1997
# Bhyb	1997
# Bhyb	2001
# Bput	2001
# Bput	2001
# Bput	2001
# Bhyb	2001
# Bhyb	2001
# Bput	2001
# Bhyb	2001
# Bhyb	2001
# Bput	2001
# Bput	2001
# Bhyb	2001
# Bhyb	2001
# Bhyb	2001
# Bhyb	2001
# Bput	2001
# Bput	2001
year = as.matrix(vegdist(metadata$Year, method = "euclidean"))

species <- read.table("Pairwise_distances_host_SNPs_ANI.txt")

# Run Mantel test to test for correlation
ani_year = mantel(ani_matrix, year, method = "spearman", na.rm = TRUE, permutations = 9999)
ani_year

ani_species = mantel(ani_matrix, species, method = "spearman", na.rm = TRUE, permutations = 9999)
ani_species

# Plotting
my_palette <- colorRampPalette(c("#8cb1d4", "white"))(n = 30)

pdf("Heatmap_ANI_BS.pdf",
  width = 10,
   height = 10,
    pointsize = 12) # Save as pdf

# Make cluster for dendogram
distance = dist(ani_matrix, method = "euclidean")
cluster = hclust(distance, method = "complete")

# Plot heatmap
heatmap.2(ani_matrix,
          main = "", # heat map title
          keysize=1.2,              # size of the colour-key
          key.title = "" ,
          key.xlab = "Average nucleotide identity",
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(13,13),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          Rowv = as.dendrogram(cluster),
          Colv= rev(as.dendrogram(cluster)),            # turn off column clustering
          cexRow = 1.2,
          cexCol = 1.2,
          RowSideColors = c(    # can be used to have coloured bars according to some feature, e.g. location or species, on the side
           rep("#ffb219", 2),
           rep("#cc0000", 1),
           rep("#ffb219", 1),
           rep("#cc0000", 1),
           rep("#ffb219", 1),
           rep("#cc0000", 2),
           rep("#ffb219", 3),
           rep("#cc0000", 1),
           rep("#ffb219", 2),
           rep("#cc0000", 3),
           rep("#ffb219", 2),
           rep("#cc0000", 1),
           rep("#ffb219", 2),
           rep("#cc0000", 2),
           rep("#ffb219", 4),
           rep("#cc0000", 2)),
          ColSideColors = c(
            rep("#788497", 3),
            rep("#afb6c2", 1),
            rep("#788497", 1),
            rep("#afb6c2", 2),
            rep("#788497", 1),
            rep("#afb6c2", 5),
            rep("#788497", 17)
          )
          )

# Adding a colour legend for the categories
# plot.new () # enable this to just plot the legend without the plot
par(lend = 2) # makes colour blocks in legend edgy

# Create a legend for the colours used in RowSideColors. Syntax: legend = c("species/location"), col = c("colour for this species/location")
legend("bottomleft",
        legend = c("Bathymodiolus hybrid - Broken Spur", "B. puteoserpentis - Broken Spur"),
        col = c("#ffb219","#cc0000"),
        lty = 2,
        lwd = 12,
        title = "Symbionts of",
        cex = 1)
legend("bottomright",
       legend = c("1997", "2001"),
       col = c("#afb6c2","#788497"),
       lty = 2,
       lwd = 12,
       title = "Sampling year",
       cex = 1)

 dev.off()               # close the pdf device