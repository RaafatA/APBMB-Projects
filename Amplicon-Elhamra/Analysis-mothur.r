#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)


sharedfile = "mothur.shared"
taxfile = "mothur.cons.taxonomy"
mapfile = "Metadata.tsv"

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

map <- read.delim(mapfile)
map <- sample_data(map)
rownames(map) <- map$group

mothur_merge <- merge_phyloseq(mothur_data, map)


colnames(tax_table(mothur_merge)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family", "Genus")

lake <- mothur_merge



sample_sum_df <- data.frame(sum = sample_sums(lake_scale))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Scale reads to even depth 
#the function
scale_reads <- function(physeq,n){
  physeq.scale <- transform_sample_counts(physeq, function(x) {round((n*x/sum(x)))})
  #otu_table(physeq.scale) = round(otu_table(physeq.scale))
  physeq.scale = prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}
lake_scale <- scale_reads(lake,min(sample_sums(lake)))


#### PCoA ####
# Ordinate
lake_pcoa <- ordinate(
  physeq = lake_scale, 
  method = "PCoA", 
  distance = "bray")
# Plot
plot_ordination(
  physeq = lake_scale,
  ordination = lake_pcoa,
  color = "depth",
  shape = "vegetation",
  title = "PCoA of Lake lake bacterial Communities"
) + 
  geom_point(aes(color = depth), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


### NMDS ###

set.seed(1)
# Ordinate
lake_nmds <- ordinate(
  physeq = lake_scale, 
  method = "NMDS", 
  distance = "bray")
plot_ordination(
  physeq = lake_scale,
  ordination = lake_nmds,
  color = "vegetation",
  shape = "depth",
  title = "NMDS of Lake lake bacterial Communities") + 
  geom_point(aes(color = vegetation), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

# Permanova

set.seed(1)
# Calculate bray curtis distance matrix
lake_bray <- phyloseq::distance(lake_scale, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(lake))
# Adonis test
adonis2(lake_bray ~ vegetation, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(lake_bray, sampledf$vegetation)
permutest(beta)



bray_cap <- phyloseq::distance(physeq = lake_scale, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = lake_scale, 
  method = "CAP",
  distance = bray_cap,
  formula = ~ pH + EC + HCO3 + Cl + SO4 + Ca + Mg + Na + K
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = lake_scale, 
  ordination = cap_ord, 
  color = "depth", 
  axes = c(1,2)
) + 
  aes(shape = vegetation) + 
  geom_point(aes(colour = depth), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5)


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

anova(cap_ord)



#######################

("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", 
 "MDS", "PCoA")

# CAP ordinate
cap_ord <- ordinate(
  physeq = lake_scale, 
  method = "CAP",
  distance = bray_cap,
  formula = ~ Cl + SO4 + Ca + Mg + Na + K + pH + HCO3 + EC
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = lake_scale, 
  ordination = cap_ord, 
  color = "depth", 
  axes = c(1,2)
) + 
  aes(shape = vegetation) + 
  geom_point(aes(colour = depth), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5)


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )


