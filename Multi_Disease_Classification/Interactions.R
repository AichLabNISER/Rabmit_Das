# Load libraries
library(dplyr)
library(readr)

# Read the files
disease_metabolite <- read_csv("D:/Downloads/Disease_Metabolite.csv")
metabolite_enzyme <- read_csv("D:/Downloads/Metabolite_Enzyme.csv")
# Clean and standardize column names
disease_metabolite <- disease_metabolite %>%
  rename(
    Disease = 1,
    Metabolites = 2,
    Interaction = 3
  )
# Merge based on the Pathway column
merged_table <- disease_metabolite %>%
  left_join(metabolite_enzyme, by = "Metabolites") %>%
  select(Disease, Enzyme, Interaction)

Enzyme_summary <- merged_table %>%
  group_by(Disease, Enzyme) %>%
  summarise(
    net_interaction = sum(Interaction),
    upregulated = sum(Interaction == 1),
    downregulated = sum(Interaction == -1),
    .groups = "drop"
  )

# Save the result
write_csv(Enzyme_summary, "D:/disease_enzyme_interaction.csv")

enzyme_pathway <- read_csv("D:/Downloads/Enzyme_Pathway.csv")
# Clean and standardize column names
disease_enzyme <- Enzyme_summary %>%
  rename(
    Disease = 1,
    Enzyme = 2,
    Interaction = 3
  )
# Merge based on the Pathway column
merged_table <- disease_enzyme %>%
  left_join(enzyme_pathway, by = "Enzyme") %>%
  select(Disease, Pathway, Interaction)

Pathway_summary <- merged_table %>%
  group_by(Disease, Pathway) %>%
  summarise(
    net_interaction = sum(Interaction),
    upregulated = sum(Interaction == 1),
    downregulated = sum(Interaction == -1),
    .groups = "drop"
  )

# Save the result
write_csv(Pathway_summary, "D:/disease_pathway_interaction.csv")

pathway_class <- read_csv("D:/Downloads/Pathway_Class.csv")
# Clean and standardize column names
disease_pathway <- Pathway_summary %>%
  rename(
    Disease = 1,
    Pathway = 2,
    Interaction = 3
  )
# Merge based on the Pathway column
merged_table <- disease_pathway %>%
  left_join(pathway_class, by = "Pathway") %>%
  select(Disease, Class, Interaction)

Class_summary <- merged_table %>%
  group_by(Disease, Class) %>%
  summarise(
    net_interaction = sum(Interaction),
    upregulated = sum(Interaction == 1),
    downregulated = sum(Interaction == -1),
    .groups = "drop"
  )

# Save the result
write_csv(Class_summary, "D:/disease_class_interaction.csv")
