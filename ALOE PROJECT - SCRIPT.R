# Preliminary Analysis Report on the Effect of Aloe arborescens Compounds on Candida albicans Morphogenesis
## MICHAEL OKANTA 


# Introduction


# Load necessary libraries
library(dplyr)
library(readr)
library(tidyverse)
library(emmeans)
library(reshape2)
library(knitr)

# Data Summary


#### Treatment Types Analyzed

# Descriptive Analysis

## Data Cleaning



# Load the data (adjust path as necessary)
aloe <- read_csv("LeeDataFungalMorphogenesis.csv")

# Rename columns for easier access (modify based on your actual column names)
colnames(aloe) <- c("Treatment", "Replicate", "Slide_Number", 
                    "No_Morphology_Change", "GT_Only", "PH_Only", 
                    "Buds_Only", "GT_PH", "GT_Buds", 
                    "PH_Buds", "Multiple_GTs")

# Convert Treatment to a factor for better handling in analysis
aloe$Treatment <- as.factor(aloe$Treatment)

# Calculate Total_Hyphae as the sum of columns where GT appears
aloe <- aloe %>%
  mutate(Total_Hyphae = GT_Only + GT_PH + GT_Buds + Multiple_GTs)

heatmap_data <- aloe %>%
  group_by(Treatment) %>%
  summarise(
    No_Morphology_Change = sum(No_Morphology_Change),
    Total_Hyphae = sum(Total_Hyphae),
    PH_Only = sum(PH_Only),
    Buds_Only = sum(Buds_Only),
    PH_Buds = sum(PH_Buds)
  ) 

heatmap_data_melt <- melt(heatmap_data, id.vars = "Treatment")

ggplot(heatmap_data_melt, aes(x = variable, y = Treatment, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Heatmap of Morphology Counts Across Treatments",
       x = "Morphology Type",
       y = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summarize the morphological responses by treatment type
morphology_summary <- aloe %>%
  group_by(Treatment) %>%
  summarise(
    No_Morphology_Change = sum(No_Morphology_Change),
    Total_Hyphae = sum(Total_Hyphae),
    PH_Only = sum(PH_Only),
    Buds_Only = sum(Buds_Only),
    PH_Buds = sum(PH_Buds)
  ) %>%
  pivot_longer(cols = -Treatment, names_to = "Morphology_Type", values_to = "Count")

# Custom color palette
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                   "#9467bd")

# Plotting the distribution of morphological responses for each treatment
ggplot(morphology_summary, aes(x = Treatment, y = Count, fill = Morphology_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Distribution of Morphological Responses Across Treatments",
       x = "Treatment",
       y = "Count of Morphological Forms") +
  theme_minimal(base_size = 14) +  # Increased base font size for readability
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),  # Rotated and adjusted text size
    plot.margin = margin(30, 30, 30, 30),  # Adds more space around the plot
    legend.position = "right"  # Positions legend on the right for clarity
  ) +
  scale_fill_manual(values = custom_colors, name = "Morphology Type") +
  theme(plot.title = element_text(size = 16, face = "bold"))


# Summary statistics for each morphology type by treatment
summary_statistics <- aloe %>%
  group_by(Treatment) %>%
  summarise(
    No_Morphology_Change_Mean = mean(No_Morphology_Change),
    No_Morphology_Change_SD = sd(No_Morphology_Change),
    No_Morphology_Change_Median = median(No_Morphology_Change),
    
    Total_Hyphae_Mean = mean(Total_Hyphae),
    Total_Hyphae_SD = sd(Total_Hyphae),
    Total_Hyphae_Median = median(Total_Hyphae),
    
    PH_Only_Mean = mean(PH_Only),
    PH_Only_SD = sd(PH_Only),
    PH_Only_Median = median(PH_Only),
    
    Buds_Only_Mean = mean(Buds_Only),
    Buds_Only_SD = sd(Buds_Only),
    Buds_Only_Median = median(Buds_Only),
    
    PH_Buds_Mean = mean(PH_Buds),
    PH_Buds_SD = sd(PH_Buds),
    PH_Buds_Median = median(PH_Buds),
    
  )

# Present`summary_statistics` as your table of summary stats
summary_statistics %>%
  kable(digits = 2,  # Round numbers to 2 decimal places
        caption = "Summary Statistics of Morphological Responses by Treatment",
        col.names = c("Treatment", 
                      "No Morphology Change Mean", "SD", "Median",
                      "Total Hyphae Mean", "SD", "Median",
                      "PH Only Mean", "SD", "Median",
                      "Buds Only Mean", "SD", "Median",
                      "PH_Buds Mean", "SD", "Median"))  # Custom column names

# Statistical Analysis Plan

  ### Model Formula:
   
## Model Diagnostics Checks 
aloe$log_offset <- log(200)

# Quasi-Poisson model with an offset
quasipoisson_model <- glm(Total_Hyphae ~ Treatment + offset(log_offset), family = quasipoisson, data = aloe)
summary(quasipoisson_model)


# Calculate the dispersion statistic for quasi-Poisson
quasi_dispersion <- sum(residuals(quasipoisson_model, type = "pearson")^2) / df.residual(quasipoisson_model)

# Pearson and deviance residuals for quasi-Poisson
pearson_resid <- residuals(quasipoisson_model, type = "pearson")
deviance_resid <- residuals(quasipoisson_model, type = "deviance")

# Residuals vs. Fitted values plot
plot(fitted(quasipoisson_model), pearson_resid, 
     xlab = "Fitted Values", ylab = "Pearson Residuals",
     main = "Residuals vs Fitted Values (Quasi-Poisson Model)")
abline(h = 0, col = "red")

# Q-Q plot for Pearson residuals
qqnorm(pearson_resid)
qqline(pearson_resid, col = "red")

# Chi-square goodness-of-fit test for quasi-Poisson model
gof_pvalue <- pchisq(sum(pearson_resid^2), df.residual(quasipoisson_model), lower.tail = FALSE)

# Calculate Cook's distance
cooksd <- cooks.distance(quasipoisson_model)

# Plot Cook's distance
plot(cooksd, type = "h", main = "Cook's Distance (Quasi-Poisson Model)", ylab = "Cook's Distance")
abline(h = 4 / length(cooksd), col = "red")

# Results & Analysis
aloe$Treatment <- factor(aloe$Treatment, 
                         levels = c("No Treatment (positive control)",
                                    "Whole Aloe arborescens Extract",
                                    "compound 'A' isolated from Aloe arborescens",
                                    "compound 'A-SA' from Sigma Aldrich",
                                    "compound 'B' isolated from Aloe arborescens",
                                    "compound 'B-SA' from Sigma Aldrich",
                                    "combination of A and B",
                                    "combination of A-SA and B-SA"), 
                         labels=c("Control", "Whole", "A", "A-SA", "B", "B-SA", "A_B", "A-SA_B-SA"))




quasipoisson_model <- glm(Total_Hyphae ~ Treatment + offset(log_offset), 
                          family = quasipoisson, data = aloe)


emms <- emmeans(quasipoisson_model, ~ Treatment)
dunnett_results <- contrast(emms, method = "trt.vs.ctrl", adjust = "dunnett")
summary(dunnett_results)


# Define the contrast for Compound A treatments (A and A-SA) vs. Compound B treatments (B and B-SA)
contrast_results <- contrast(emms, 
                             list(A_vs_B = c(0, 0, 0.5, 0.5, -0.5, -0.5, 0, 0)),
                             adjust = "none")

# Display the contrast results
summary(contrast_results)


# Define contrasts to compare Whole Aloe extract with the six other treatments (excluding Control)
# Specify contrasts comparing "Whole" with each other treatment (excluding "Control")
comparison_results <- contrast(emms, method = "trt.vs.ctrl", ref = "Whole", 
                               exclude = "Control", adjust = "none")

# Display the results
summary(comparison_results)

# Conclusion
