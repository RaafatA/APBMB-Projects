library(vegan)
library(tidyverse)
library(microeco)
library(magrittr)

setwd("/home/raafat/Documents/16s Amplicon/CCA")
count_matrix <- read.csv("otu_table.csv", header = TRUE, row.names = 1)
total_abundance <- colSums(count_matrix)
relative_abundance <- count_matrix/total_abundance


data("env_data_16S")

OTO_TABLE = read.csv("otu_table.csv", row.names = 1)
tax_table = read.csv("taxonomy.csv", row.names = 1)
sample_info = read.csv("sample_info.csv", row.names = 1)
env_data = read.csv("env_info__scaled.csv", row.names = 1)
OTO_TABLE %<>% cal
tax_table %<>% tidy_taxonomy

#
#create a microtable object
#dataset <- microtable$new(otu_table = OTO_TABLE)
# generally add the metadata
#dataset <- microtable$new(otu_table = OTO_TABLE, sample_table = sample_info,tax_table = tax_table)
#dataset
#

#create a microtable object
object1 <- microtable$new(sample_table = sample_info, otu_table = OTO_TABLE, tax_table = tax_table)
object1
env_data[,"Na"]
# add_data is used to add the environmental data
t1 <- trans_env$new(dataset = object1, add_data = env_data[,1:4])

new_test <- clone(object1)
new_test$sample_table <- data.frame(new_test$sample_table, env_data[rownames(new_test$sample_table), ])
# now new_test$sample_table has the whole data
new_test
t1 <- trans_env$new(dataset = new_test, env_cols = 6:14, standardize = TRUE)








# use Wilcoxon Rank Sum Test as an example
t1$cal_diff(group = "site", method = "wilcox")
t1$data_env
head(t1$res_diff)



# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = relative_abundance, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
dataset



t1$cal_diff(method = "anova", group = "site")
# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}
plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 3))


t1$cal_diff(group = "depth", by_group = "site", method = "anova")
t1$plot_diff(measure = "pH", add_sig_text_size = 5)


t1$cal_autocor()
t1$cal_autocor(group = "depth")


# use bray-curtis distance for dbRDA
t1$cal_ordination(method = "dbRDA", use_measure = "bray")
# t1$res_rda is the result list stored in the object
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.2)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "site", plot_shape = "depth")


data(env_data_16S)
# Caclulating the CCA
t1$cal_ordination(method = "CCA", taxa_level = "Phylum")
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE)
t1$cal_ordination_anova()
t1$cal_ordination_envfit()
t1$plot_ordination(plot_color = "site", plot_shape = "depth")
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
head(t1$res_mantel)
t1$trans_ordination(adjust_arrow_length = TRUE)
g1 <- t1$plot_ordination(plot_color = "site", plot_shape = "depth")

ggplot2::ggsave("RDA.svg", g1, width = 8, height = 6.5)
# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_ordination(show_taxa = 5, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 1, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
t1$plot_ordination(plot_color = "site")


t1$cal_autocor(group = "site")



#-----------------------------------------------------------------

# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))

tabun <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
tabun$plot_bar(others_color = "grey70", facet = "depth", xtext_keep = FALSE, legend_text_italic = FALSE)
tabun$plot_box(group = "site", xtext_angle = 30)




new_test <- clone(object1)
new_test$sample_table <- data.frame(new_test$sample_table, env_data[,1:10])
new_test$sample_table
# now new_test$sample_table has the whole data
new_test
t1 <- trans_env$new(dataset = new_test, env_cols = 6:15,standardize = TRUE)
# use bray-curtis distance for dbRDA
t1$cal_ordination(method = "dbRDA", use_measure = "bray")
# t1$res_rda is the result list stored in the object
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = .5)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "site", plot_shape ="depth")
data(env_data_16S)
#################################################################################
#############################################################################3##
new_test <- clone(object1)
new_test$sample_table <- data.frame(new_test$sample_table, env_data[rownames(new_test$sample_table), ])
# now new_test$sample_table has the whole data
new_test
# now new_test$sample_table has the whole data
