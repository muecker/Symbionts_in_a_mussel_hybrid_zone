# Plot heatmap of Bathymodiolus SOX symbiont FST values and test correlation with sampling year & host species 2020-01-15
# Merle Ücker

# Load libraries
library(gplots)
library(RColorBrewer)
library(vegan)

# Load data
fst <- read.table("Pairwise_mean_FST_SOX_symbionts_Broken_Spur.csv", sep = "\t", header = TRUE, check.names = FALSE)

# Set row names to bin names
rownames <- fst[,1]
row.names(fst) <- rownames
fst[,1] <- NULL
fst_matrix <- as.matrix(fst)

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

species <- read.table("Pairwise_distances_host_SNPs_FST.txt")

# Run Mantel test to test for correlation
fst_year = mantel(fst_matrix, year, method = "spearman", na.rm = TRUE, permutations = 9999)
fst_year

fst_species = mantel(fst_matrix, species, method = "spearman", na.rm = TRUE, permutations = 9999)
fst_species

# Plotting
my_palette <- colorRampPalette(c("white", "gray88", "#00bcd4"))(n = 30)

# Save as pdf image
pdf("Heatmap_FST_BS_nodensity.pdf",
  width = 10,
   height = 10,
    pointsize = 12)

# Make cluster for dendogram
distance = dist(fst_matrix, method = "euclidean")
cluster = hclust(distance, method = "complete")

# Plot heatmap
heatmap.2(fst_matrix,
          main = "", # heat map title
          keysize=1.2,
          key.title = "" ,
          key.xlab = expression('F'[ST]),
          key.ylab = "Density",
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          Rowv = as.dendrogram(cluster),
          Colv= rev(as.dendrogram(cluster)),            # turn off column clustering
          cexRow = 1.2,
          cexCol = 1.2,
          RowSideColors = c(    # can be used to have coloured bars according to some feature, e.g. location or species, on the side
            rep("#ffb219", 2),
            rep("#cc0000", 2),
            rep("#ffb219", 1),
            rep("#cc0000", 1),
            rep("#ffb219", 4),
            rep("#cc0000", 1),
            rep("#ffb219", 2),
            rep("#cc0000", 3),
            rep("#ffb219", 2),
            rep("#cc0000", 1),
            rep("#ffb219", 3),
            rep("#cc0000", 2),
            rep("#ffb219", 4),
            rep("#cc0000", 2)),
          ColSideColors = c(
            rep("#788497", 1),
            rep("#afb6c2", 1),
            rep("#788497", 1),
            rep("#afb6c2", 9),
            rep("#788497", 18)))

# Adding a colour legend for the categories
# plot.new () # enable this to just plot the legend without the plot
par(lend = 2) # makes colour blocks in legend edgy

# Create a legend for the colours used in RowSideColors. Syntax: legend = c("species/location"), col = c("colour for this species/location")
legend("bottomleft",
        legend = c("Bathymodiolus hybrid", "B. puteoserpentis"),
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

 dev.off()               # close the PDF device