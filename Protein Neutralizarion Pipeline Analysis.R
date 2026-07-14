#--------------------------------------------
          Demonstration version 
#--------------------------------------------

### HIV Protein Neutralization Analysis #####

### Date Created: 01-14-2026 ####

### Author: D'Atra Hill####

### Description:Analysis pipeline for MTT cell viability and β-gal infectivity 
#assays using R, tidyverse, and ggplot2.

### © 2026 D'Atra Hill. All rights reserved.

#-------------------------------------------
##Upload Data
# Load libraries

library(ggplot2)
library(dplyr)
library(tidyr)

                                        # --- Input Raw MTT data ---
MTT_data <- data.frame(
  Blank = c(0.518, 0.462, 0.497, 0.438, 0.524, 0.483),      
  Cells = c(0.887, 1.108, 0.791, 0.953, 0.846, 0.999),      
  control = c(0.614, 0.553, 0.596, 0.571, 0.542, 0.629), 
  Dilution_1 = c(0.479, 0.457, 0.454, 0.445, 0.420, 0.447),
  Dilution_2 = c(0.531, 0.489, 0.445, 0.479, 0.458, 0.502),
  Dilution_3 = c(0.944, 1.017, 1.184, 0.865, 1.126, 1.271)
)

# --- Average blank ---
avg_blank <- mean(MTT_data$Blank)
# --- Correct samples by subtracting blank ---
MTT_corrected <- MTT_data %>%
  mutate(across(Cells:Dilution_3, ~ . - avg_blank)) %>%
  summarise(across(Cells:Dilution_3, mean))
# --- Prepare MTT plot data ---
MTT_plot <- data.frame(
  Sample = c("Blank", "Cells", "Control", "Dilution_1", "Dilution_2", "Dilution_3"),
  MTT_OD = c(avg_blank, as.numeric(MTT_corrected[c("Cells","Control","Dilution_1","Dilution_2","Dilution_3")]))
)

# --- Ensure factor levels are ordered ---
MTT_plot$Sample <- factor(MTT_plot$Sample, levels = MTT_plot$Sample)

                                    # --- MTT plot ---

# --- MTT % viability plot ---
ggplot(MTT_percent_plot, aes(x = Sample, y = Percent_Viability, fill = Sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percent_Viability,1), "%")), vjust=-0.5) +
  labs(title = "MTT Assay (% Cell Viability)",
       x = "Sample",
       y = "% Cell Viability") +
  scale_fill_manual(values = c(
    "Blank"="#d9d9d9",
    "Cells"="darkgray",
    "Control"="goldenrod1",
    "Dilution_1"="lightpink",
    "Dilution_2"="deeppink3",
    "Dilution_3"="darkmagenta")) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Combine all samples including blank
MTT_corrected_full <- c(avg_blank, as.numeric(MTT_corrected[c("Cells","Control","Dilution_1","Dilution_2","Dilution_3")]))
# Calculate % viability relative to TZM-bl (second element)
percent_viability <- (MTT_corrected_full / MTT_corrected_full[2]) * 100
# Prepare plot data
MTT_percent_plot <- data.frame(
  Sample = c("Blank", "Cells", "Control", "Dilution_1", "Dilution_2", "Dilution_3"),
  Percent_Viability = percent_viability
)

                                  # Order factor for plotting
MTT_percent_plot$Sample <- factor(MTT_percent_plot$Sample, levels = MTT_percent_plot$Sample)
# --- MTT % viability plot ---
ggplot(MTT_percent_plot, aes(x = Sample, y = Percent_Viability, fill = Sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percent_Viability,1), "%")), vjust=-0.5) +
  labs(title = "MTT Assay (% Cell Viability)",
       x = "Sample",
       y = "% Cell Viability") +
  scale_fill_manual(values = c(
    "Blank"="#d9d9d9",
    "Cells"="darkgray",
    "Control"="goldenrod1",
    "Dilution_1"="lightpink",
    "Dilution_2"="deeppink3",
    "Dilution_3"="darkmagenta")) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
                                      #---Raw B-Gal data---

beta_gal_data <- data.frame(
  Cells= c(0.518, 0.462, 0.497, 0.438, 0.524, 0.483),
  Control = c(0.887, 1.108, 0.791, 0.953, 0.846, 0.999),
  Dilution_1 = c(0.614, 0.553, 0.596, 0.571, 0.542, 0.629),
  Dilution_2 = c(0.531, 0.489, 0.445, 0.479, 0.458, 0.502),
  Dilution_3 = c(0.944, 1.017, 1.184, 0.865, 1.126, 1.271)
)

#----- Correct Beta-Gal for baseline (cells_only)-------
beta_gal_corrected <- beta_gal_data %>%
  mutate(across(Protein_control:Dilution_3, ~. - Cells_only)) %>%
  summarise(across(Protein_control:Dilution_3, mean))

# --- Calculate % viral inhibition ---
beta_gal_inhibition <- 1 - (beta_gal_corrected[, c("Dilution_2","Dilution_2","Dilution_3")] /
                              beta_gal_corrected$Protein_control)
beta_gal_inhibition <- as.numeric(beta_gal_inhibition) * 100

# --- Prepare β-gal plot data in order with short labels ---
beta_plot <- data.frame(
  Sample_short = c("Cells","Control","Dilution_2","Dilution_2","Dilution_3"),
  Viral_Inhibition = c(0, 0, beta_gal_inhibition)
)

# --- Ensure factor order ---
beta_plot$Sample_short <- factor(beta_plot$Sample_short, levels = beta_plot$Sample_short)
# --- β-gal plot ---
ggplot(beta_plot, aes(x = Sample_short, y = Viral_Inhibition, fill = Sample_short)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Viral_Inhibition,1), "%")), vjust=-0.5) +
  labs(title = "β-gal Assay (% Viral Inhibition) with Controls",
       x = "Sample",
       y = "% Viral Inhibition") +
  scale_fill_manual(values = c("Cells"="darkgrey","Control"="goldenrod1",
                               "Dilution_3"="darkmagenta","Dilution_1"="lightpink","Dilution_2"="deeppink3")) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#---Raw B-Gal data---
beta_gal_data <- data.frame(
  Cells= c(0.518, 0.462, 0.497, 0.438, 0.524, 0.483),
  Control = c(0.887, 1.108, 0.791, 0.953, 0.846, 0.999),
  Dilution_1 = c(0.614, 0.553, 0.596, 0.571, 0.542, 0.629),
  Dilution_2 = c(0.531, 0.489, 0.445, 0.479, 0.458, 0.502),
  Dilution_3 = c(0.944, 1.017, 1.184, 0.865, 1.126, 1.271)
)

# --- Calculate % infectivity relative to Protein control ---
beta_gal_corrected_means <- beta_gal_data %>%
  mutate(across(Protein_control:Dilution_3, ~ . - Cells_only)) %>%
  summarise(across(Protein_control:Dilution_3, mean))

# Infectivity of each treatment as % of Protein control
beta_gal_percent_infectivity <- (beta_gal_corrected_means / beta_gal_corrected_means$Protein_control) * 100
# Prepare data for plotting
beta_plot <- data.frame(
  Sample = c("Cells", "Control", "Dilution_1", "Dilution_2", "Dilution_3"),
  Percent_Infectivity = c(0, 100,
                          beta_gal_percent_infectivity$Dilution_1,
                          beta_gal_percent_infectivity$Dilution_2,
                          beta_gal_percent_infectivity$Dilution_3)
)

# Ensure correct order
beta_plot$Sample <- factor(beta_plot$Sample, levels = beta_plot$Sample)
# --- Plot β-gal infectivity ---
ggplot(beta_plot, aes(x = Sample, y = Percent_Infectivity, fill = Sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percent_Infectivity,1), "%")), vjust=-0.5) +
  labs(title = "β-gal Assay (% Infectivity Relative to Protein Control)",
       x = "Sample",
       y = "% Infectivity") +
  scale_fill_manual(values = c(
    "Cells"="darkgrey",
    "Control"="goldenrod1",
    "Dilution_1"="lightpink",
    "Dilution_2"="deeppink3",
    "Dilution_3"="darkmagenta"
  )) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
#---Raw B-Gal data---
beta_gal_data <- data.frame(
  Cells= c(0.518, 0.462, 0.497, 0.438, 0.524, 0.483),
  Control = c(0.887, 1.108, 0.791, 0.953, 0.846, 0.999),
  Dilution_1 = c(0.614, 0.553, 0.596, 0.571, 0.542, 0.629),
  Dilution_2 = c(0.531, 0.489, 0.445, 0.479, 0.458, 0.502),
  Dilution_3 = c(0.944, 1.017, 1.184, 0.865, 1.126, 1.271)
)
### Date Created: 02-19-2026 ####
### Creator Name: D'Atra Hill####
### Description:Corrections and added addtions based on data
### © 2026 D'Atra Hill. All rights reserved.

# --- Calculate % infectivity relative to Protein control ---
beta_gal_corrected_means <- beta_gal_data %>%
  mutate(across(Protein_control:Dilution_3, ~ . - Cells_only)) %>%
  summarise(across(Protein_control:Dilution_3, mean))
# Infectivity of each treatment as % of Protein control
beta_gal_percent_infectivity <- (beta_gal_corrected_means / beta_gal_corrected_means$Protein_control) * 100
# Prepare data for plotting
beta_plot <- data.frame(
  Sample = c("Cells", "Control", "Dilution_1", "Dilution_3", "Dilution_3"),
  Percent_Infectivity = c(0, 100,
                          beta_gal_percent_infectivity$Dilution_1,
                          beta_gal_percent_infectivity$Dilution_2,
                          beta_gal_percent_infectivity$Dilution_3)
)

# Ensure correct order of the data and the code 
beta_plot$Sample <- factor(beta_plot$Sample, levels = beta_plot$Sample)
# --- Plot β-gal infectivity ---
ggplot(beta_plot, aes(x = Sample, y = Percent_Infectivity, fill = Sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percent_Infectivity,1), "%")), vjust=-0.5) +
  labs(title = "β-gal Assay (% Infectivity Relative to Protein Control)",
       x = "Sample",
       y = "% Infectivity") +
  scale_fill_manual(values = c(
    "Cells"="darkgrey",
    "Control"="goldenrod1",
    "Dilution_1"="lightpink",
    "Dilution_2"="deeppink3",
    "Dilution_3"="darkmagenta"
  )) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
#---Raw B-Gal data---
beta_gal_data <- data.frame(
  Cells= c(0.518, 0.462, 0.497, 0.438, 0.524, 0.483),
  Control = c(0.887, 1.108, 0.791, 0.953, 0.846, 0.999),
  Dilution_1 = c(0.614, 0.553, 0.596, 0.571, 0.542, 0.629),
  Dilution_2 = c(0.531, 0.489, 0.445, 0.479, 0.458, 0.502),
  Dilution_3 = c(0.944, 1.017, 1.184, 0.865, 1.126, 1.271)
)
                    # --- Calculate % infectivity relative to control ---
beta_gal_corrected_means <- beta_gal_data %>%
  mutate(across(Protein_control:Dilution_3, ~ . - Cells_only)) %>%
  summarise(across(Protein_control:Dilution_3, mean))

                  # Infectivity of each treatment as % of control
beta_gal_percent_infectivity <- (beta_gal_corrected_means / beta_gal_corrected_means$Protein_control) * 100
# Prepare data for plotting
beta_plot <- data.frame(
  Sample = c("Cells", "Control", "Dilution_1", "Dilution_2", "Dilution_3"),
  Percent_Infectivity = c(0, 100,
                          beta_gal_percent_infectivity$Dilution_1,
                          beta_gal_percent_infectivity$Dilution_2,
                          beta_gal_percent_infectivity$Dilution_3)
)
# Ensure correct order
beta_plot$Sample <- factor(beta_plot$Sample, levels = beta_plot$Sample)

                          # --- Plot β-gal infectivity ---
ggplot(beta_plot, aes(x = Sample, y = Percent_Infectivity, fill = Sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percent_Infectivity,1), "%")), vjust=-0.5) +
  labs(title = "β-gal Assay (% Infectivity Relative to Protein Control)",
       x = "Sample",
       y = "% Infectivity") +
  scale_fill_manual(values = c(
    "Cells"="darkgrey",
    "Control"="goldenrod1",
    "Dilution_1"="lightpink",
    "Dilution_2"="deeppink3",
    "Dilution_3"="darkmagenta"
  )) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# --- Raw MTT data (replace with your measured OD values) ---
mtt_data <- tribble(
  ~Blank, ~Cells, ~Protein_control, ~Dilution_1, ~Dilution_2, ~Dilution_3,
  0.518, 0.887, 0.614, 0.531, 0.944,
  0.462, 1.108, 0.553, 0.489, 1.017,
  0.497, 0.791, 0.596, 0.445, 1.184,
  0.438, 0.953, 0.571, 0.479, 0.865,
  0.524, 0.846, 0.542, 0.458, 1.126,
  0.483, 0.999, 0.629, 0.502, 1.271
)
# --- Calculate % viability ---
mtt_means <- mtt_data %>% summarise(across(everything(), mean))
blank_mean <- mtt_means$Blank
cells_mean <- mtt_means$Cells
mtt_percent <- (mtt_means - blank_mean) / (cells_mean - blank_mean) * 100
mtt_percent <- mtt_percent %>% select(-Blank)
# --- Convert to long format for plotting ---
mtt_long <- mtt_percent %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "Percent_Viability")
# --- Order the conditions ---
mtt_long$Condition <- factor(mtt_long$Condition,
                             levels = c("Cells", "Dilution_1", "Dilution_2", "Dilution_3"))
# --- Plot ---
ggplot(mtt_long, aes(x = Condition, y = Percent_Viability, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(Percent_Viability,1)), vjust=-0.5) +
  scale_fill_manual(values = c("Cells"="darkgray",
                               "Dilution_1"="lightpink",
                               "Dilution_2"="deeppink3",
                               "Dilution_3"="darkmagenta")) +
  labs(title="MTT Assay: % Cell Viability", y="Cell Viability (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- Raw MTT data ---
#This will provide violin plots to visualize the date fully with errors of distirbution
mtt_data <- tribble(
  ~Blank, ~Cells, ~Protein_control, ~Dilution_1, ~Dilution_2, ~Dilution_3,
  0.518, 0.887, 0.614, 0.531, 0.944,
  0.462, 1.108, 0.553, 0.489, 1.017,
  0.497, 0.791, 0.596, 0.445, 1.184,
  0.438, 0.953, 0.571, 0.479, 0.865,
  0.524, 0.846, 0.542, 0.458, 1.126,
  0.483, 0.999, 0.629, 0.502, 1.271
)
# --- Normalize for each replicate ---
mtt_normalized <- mtt_data %>%
  mutate(across(-Blank, ~ (.x - Blank))) %>%  # subtract blank for each sample
  mutate(across(-Cells, ~ .x / (Cells - Blank) * 100)) # normalize to cells mean = 100%
# --- Summarize by condition ---
mtt_summary <- mtt_normalized %>%
  select(-Blank) %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "Percent_Viability") %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Percent_Viability, na.rm = TRUE),
    SD = sd(Percent_Viability, na.rm = TRUE)
  )
# --- Order conditions ---
mtt_summary$Condition <- factor(
  mtt_summary$Condition,
  levels = c("Cells", "Protein_control", "Dilution_1", "Dilution_2", "Dilution_3")
)
# --- Plot ---
ggplot(mtt_summary, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  geom_text(aes(label = round(Mean, 1)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(
    "Cells" = "darkgray",
    "Protein_control" = "goldenrod1",
    "Dilution_1" = "lightpink",
    "Dilution_2" = "deeppink3",
    "Dilution_3" = "darkmagenta"
  )) +
  labs(
    title = "MTT Assay: Percent Cell Viability",
    y = "Cell Viability (%)",
    x = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# --- Raw MTT data ------
mtt_data <- tribble(
  ~Blank, ~Cells, ~Protein_control, ~Dilution_1, ~Dilution_2, ~Dilution_3,
  0.518, 0.887, 0.614, 0.531, 0.944,
  0.462, 1.108, 0.553, 0.489, 1.017,
  0.497, 0.791, 0.596, 0.445, 1.184,
  0.438, 0.953, 0.571, 0.479, 0.865,
  0.524, 0.846, 0.542, 0.458, 1.126,
  0.483, 0.999, 0.629, 0.502, 1.271
)
# --- Calculate % viability ---
mtt_means <- mtt_data %>% summarise(across(everything(), mean))
blank_mean <- mtt_means$Blank
cells_mean <- mtt_means$Cells
mtt_percent <- (mtt_means - blank_mean) / (cells_mean - blank_mean) * 100
mtt_percent <- mtt_percent %>% select(-Blank)
# --- Convert to long format for plotting ---
mtt_long <- mtt_percent %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "Percent_Viability")
# --- Order the conditions ---
mtt_long$Condition <- factor(mtt_long$Condition,
                             levels = c("Cells", "Dilution_1", "Dilution_2", "Dilution_3"))
# --- Plot ---
ggplot(mtt_long, aes(x = Condition, y = Percent_Viability, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(Percent_Viability,1)), vjust=-0.5) +
  scale_fill_manual(values = c("Cells"="darkgray",
                               "Dilution_1"="lightpink",
                               "Dilution_2"="deeppink3",
                               "Dilution_3"="darkmagenta")) +
  labs(title="MTT Assay: % Cell Viability", y="Cell Viability (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- Raw MTT data 
mtt_data <- tribble(
  ~Blank, ~Cells, ~Control, ~Dilution_1, ~Dilution_2, ~Dilution_3,
  0.518, 0.887, 0.614, 0.531, 0.944,
  0.462, 1.108, 0.553, 0.489, 1.017,
  0.497, 0.791, 0.596, 0.445, 1.184,
  0.438, 0.953, 0.571, 0.479, 0.865,
  0.524, 0.846, 0.542, 0.458, 1.126,
  0.483, 0.999, 0.629, 0.502, 1.271
)
# --- Calculate % viability ---
mtt_means <- mtt_data %>% summarise(across(everything(), mean))
blank_mean <- mtt_means$Blank
cells_mean <- mtt_means$Cells
mtt_percent <- (mtt_means - blank_mean) / (cells_mean - blank_mean) * 100
mtt_percent <- mtt_percent %>% select(-Blank)
# --- Convert to long format for plotting ---
mtt_long <- mtt_percent %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "Percent_Viability")
# --- Order the conditions ---
mtt_long$Condition <- factor(mtt_long$Condition,
                             levels = c("Cells","Control" ,"Dilution_1", "Dilution_2", "Dilution_3"))
# --- Plot ---
ggplot(mtt_long, aes(x = Condition, y = Percent_Viability, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(Percent_Viability,1)), vjust=-0.5) +
  scale_fill_manual(values = c("Cells"="darkgray",
                               "Control"="goldenrod1",
                               "Dilution_1"="lightpink",
                               "Dilution_2"="deeppink3",
                               "Dilution_3"="darkmagenta")) +
  labs(title="MTT Assay: % Cell Viability", y="Cell Viability (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


