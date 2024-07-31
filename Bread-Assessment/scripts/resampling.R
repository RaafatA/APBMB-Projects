installed.packages(c("moderndive","infer"))
library(tidyverse)
library(moderndive)
library(infer)
library(ggpubr)
setwd("/home/raafat/Documents/Bread_Wheat")
data = read.csv("/home/raafat/Documents/Bread_Wheat/data/Planting Date 2 GH Agronomic Traits - Sheet1.csv")
library(infer)

# Specify the data and variable of interest
specify_data <- specify(data, response = LL)

# Generate bootstrap samples
bootstrap_samples <- generate(specify_data, reps = 1000, type = "bootstrap") %>% calculate(stat = "mean")

na_indices <- which(is.na(data$LA))

# For each NA value
for(i in na_indices) {
  
  data$LA[i] <- sample(bootstrap_samples$stat, 1)
}

write.csv(data, file = "/home/raafat/Documents/Bread_Wheat/data/Planting Date 2 GH Agronomic Traits - Sheet1_resampled",quote = FALSE, sep = ",")
write.csv()




# Analysis ########################3
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
# Load necessary library
library(xtable)
library(gridExtra)  # for arranging plots

# Perform ANOVA for each trait and print a publication-ready table
for (trait in traits) {
  aov_result <- aov(as.formula(paste(trait, "~ genotype + treatment")), data = data)
  anova_table <- Anova(aov_result, type="II")
  
  # Convert the ANOVA table to LaTeX
  latex_table <- xtable(anova_table, caption = paste("ANOVA Table for", trait))
  print(latex_table, type = "latex", file = paste0(trait, "_anova_table.tex"))
}
gather(Treatment)


# Calculate standard error for each group
data_se <- Data %>%
  gather(treatment,PH) %>%
  group_by(Treatment, PH) %>%
  summarise(value_mean = mean(value), se = sd(value) / sqrt(n()))

# Plot each trait with error bars
ggplot(data_se, aes(x=interaction(genotype, treatment), y=value_mean, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) +
  geom_errorbar(aes(ymin=value_mean-se, ymax=value_mean+se), width=0.2, position=position_dodge(0.7)) +
  facet_wrap(~trait, scales="free") +
  labs(x="Genotype and Treatment", y="Value") +
  theme_minimal()

# Reshape the data
df_long <- data %>%
  gather(key = "Trait", value = "Value", PH, LL, LW, LA, SL)

# Calculate means and standard errors
summary <- df_long %>%
  group_by(genotype, treatment, Trait) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n())
  )
plots <- list()

# Create a bar plot for each trait
for(trait in unique(df_long$Trait)) {
  p <- ggplot(summary[summary$Trait == trait,], aes(x = treatment, y = mean, fill = summary$genotype[1])) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
    scale_fill_brewer(palette = "Set2") +  # Use a different color palette
    labs(title = paste("Trait:", trait), x = "Treatment", y = "Value") +
    theme_minimal()
  plots[[trait]] <- p
  print(p)
}
do.call(grid.arrange, plots)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)

# Calculate means, standard errors, and ANOVA for significance
summary <- df_long %>%
  group_by(genotype, treatment, Trait) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  group_by(Trait) %>%
  do(tidy(aov(mean ~ treatment, data = .)))

# Filter out only significant comparisons
significant_comparisons <- summary %>%
  filter(p.value < 0.05)

plots <- list()

# Create a combined bar plot for each trait
for(trait in unique(df_long$Trait)) {
  # Filter the significant comparisons for the current trait
  current_significant <- significant_comparisons %>%
    filter(Trait == trait)
  
  # Create the plot
  p <- ggplot(df_long[df_long$Trait == trait,], aes(x = treatment, y = mean, fill = genotype)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste("Trait:", trait), x = "Treatment", y = "Value") +
    theme_minimal() +
    # Add significance annotations
    geom_text(data = current_significant, aes(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", "*")),
                                              y = max(df_long$mean) * 1.1), # Adjust y position for annotations
              position = position_dodge(width = 0.75), hjust = 0.5, vjust = 0)
  
  # Save the plot in a list
  plots[[trait]] <- p
  
  # Print the plot
  print(p)
}



# Load necessary libraries
install.packages(c("dplyr", "car", "purrr"))
library(dplyr)
library(car)
library(purrr)

# Read the data from a CSV file
data <- read.csv("data.csv")

# Function to perform ANOVA on each trait
perform_anova <- function(trait) {
  formula <- as.formula(paste(trait, "~ genotype * treatment"))
  model <- aov(formula, data = data)
  summary(model)[[1]][, c("Df", "F value", "Pr(>F)")]
}

# Apply the function to each trait
results <- map(traits, perform_anova)

# Convert the list to a data frame and transpose it
results_df <- as.data.frame(do.call(cbind, results))
colnames(results_df) <- traits
results_df <- t(results_df)

# Print the results
print(results_df)










