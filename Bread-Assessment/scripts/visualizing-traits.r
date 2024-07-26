





# Load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(agricolae)
library(tidyr) 
library(purrr) 

setwd("/home/raafat/Documents/Bread_Wheat/Final_Data")
df = read.csv("Biochemical-Bread-openfield-2023.csv")

# Melt the data frame
df_melt <- df %>% 
  melt(id.vars = c("Genotype", "Treatment", "Season"), measure.vars = c("APX"))

# Calculate summary statistics
summary_df <- df_melt %>%
  group_by(Season, Cultivar, Treatment) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()))

# Merge summary statistics back with the original melted data for plotting
df_plot <- left_join(df_melt, summary_df, by = c("Season", "Cultivar", "Treatment"))

# Perform LSD test
df_lsd <- df %>% 
  group_by(Season, Cultivar) %>%
  nest() %>%
  mutate(lsd_results = map(data, ~ {
    model <- aov(FD ~ Treatment, data = .x)
    lsd <- LSD.test(model, "Treatment", group = TRUE)$groups
    lsd$Treatment <- rownames(lsd)
    return(lsd)
  })) %>%
  unnest(cols = lsd_results)

# Rename the columns for clarity
df_lsd <- df_lsd %>%
  rename(lsd_group = groups)

# Merge LSD results with the plot data
df_plot <- df_plot %>%
  left_join(df_lsd, by = c("Season", "Cultivar", "Treatment"))

# Custom colors
custom_colors <- c("T2" = "grey30", "T1" = "grey60", "Control" = "grey90")

# Function to create a plot
create_plot <- function(title) {
  ggplot(df_plot, aes(x = Cultivar, y = mean_value, fill = Treatment)) + 
    geom_col(position = position_dodge(0.8), width = 0.8, color = "black", size=0.7) + # Increase dodge width for spacing
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                  width = 0.2, position = position_dodge(0.8)) +  # Adjust dodge width for error bars
    facet_wrap(~ Season, scales = "free_x", nrow = 1, strip.position = "top") +  # Split by Season, adjust layout
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    labs(title = title, 
         x = "Cultivar", y = "Flowering Dates") + 
    theme_classic() +  # Change to classic theme to remove the grid
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Times New Roman", face = "bold", size = 12),  # Font style for x-axis labels
      axis.text.y = element_text(family = "Times New Roman", face = "bold"),  # Font style for y-axis labels
      axis.title.x = element_text(family = "Times New Roman", face = "bold"),  # Font style for x-axis title
      axis.title.y = element_text(family = "Times New Roman", face = "bold"),  # Font style for y-axis title
      plot.title = element_text(family = "Times New Roman", face = "bold"),  # Font style for plot title
      strip.background = element_blank(),  # Remove facet background
      strip.text = element_text(size = 14, family = "Times New Roman", face = "bold"),  # Font style for facet labels
      panel.spacing = unit(2, "lines"),  # Increase space between facets
      legend.text = element_text(family = "Times New Roman", face = "bold"),  # Font style for legend text
      legend.title = element_text(family = "Times New Roman", face = "bold")  # Font style for legend title
    ) +  
    geom_text(aes(label = lsd_group), 
              position = position_dodge(0.8), 
              vjust = -1.4)  # Add LSD group annotations
}

# Create multiple plots
plot1 <- create_plot("(A)")
plot2 <- create_plot("(B)")
plot3 <- create_plot("(C)")
plot4 <- create_plot("(D)")
plot5 <- create_plot("(E)")
plot6 <- create_plot("(F)")
# Arrange plots in a grid
grid.arrange(plot1, plot2,ncol = 2)








##########################################################################################################3
# Melt the data frame
df_melt <- df %>% 
  melt(id.vars = c("Genotype", "Treatment"), measure.vars = c("POD"))

# Calculate summary statistics
summary_df <- df_melt %>%
  group_by(Genotype, Treatment) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()))

# Merge summary statistics back with the original melted data for plotting
df_plot <- left_join(df_melt, summary_df, by = c("Genotype", "Treatment"))

# Perform LSD test
df_lsd <- df %>% 
  group_by(Genotype) %>%
  nest() %>%
  mutate(lsd_results = map(data, ~ {
    model <- aov(POD ~ Treatment, data = .x)
    lsd <- LSD.test(model, "Treatment", group = TRUE)$groups
    lsd$Treatment <- rownames(lsd)
    return(lsd)
  })) %>%
  unnest(cols = lsd_results)

# Rename the columns for clarity
df_lsd <- df_lsd %>%
  rename(lsd_group = groups)

# Merge LSD results with the plot data
df_plot <- df_plot %>%
  left_join(df_lsd, by = c("Genotype", "Treatment"))

# Custom colors
custom_colors <- c("T2" = "grey30", "T1" = "grey60", "C" = "grey90")

# Function to create a plot
create_plot <- function(title) {
  ggplot(df_plot, aes(x = Genotype, y = mean_value, fill = Treatment)) + 
    geom_col(position = position_dodge(0.8), width = 0.8, color = "black", size=0.7) + # Increase dodge width for spacing
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                  width = 0.2, position = position_dodge(0.8)) +  # Adjust dodge width for error bars
    #facet_wrap(~ Season, scales = "free_x", nrow = 1, strip.position = "top") +  # Split by Season, adjust layout
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    labs(title = title, 
         x = "Genotype", y = "Peroxidase ( mM/gFW )") + 
    theme_classic() +  # Change to classic theme to remove the grid
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Times New Roman", face = "bold", size = 12),  # Font style for x-axis labels
      axis.text.y = element_text(family = "Times New Roman", face = "bold"),  # Font style for y-axis labels
      axis.title.x = element_text(family = "Times New Roman", face = "bold"),  # Font style for x-axis title
      axis.title.y = element_text(family = "Times New Roman", face = "bold"),  # Font style for y-axis title
      plot.title = element_text(family = "Times New Roman", face = "bold"),  # Font style for plot title
      strip.background = element_blank(),  # Remove facet background
      strip.text = element_text(size = 14, family = "Times New Roman", face = "bold"),  # Font style for facet labels
      panel.spacing = unit(2, "lines"),  # Increase space between facets
      legend.text = element_text(family = "Times New Roman", face = "bold"),  # Font style for legend text
      legend.title = element_text(family = "Times New Roman", face = "bold")  # Font style for legend title
    ) +  
    geom_text(aes(label = lsd_group), 
              position = position_dodge(0.8), 
              vjust = -1.4)  # Add LSD group annotations
}

# Create multiple plots
plot1 <- create_plot("(A)")
plot2 <- create_plot("(B)")
plot3 <- create_plot("(C)")
plot4 <- create_plot("(D)")
plot5 <- create_plot("(E)")
plot6 <- create_plot("(F)")
# Arrange plots in a grid
grid.arrange(plot1, plot2,plot3, plot4,ncol = 2)

