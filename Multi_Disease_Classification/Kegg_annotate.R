# -----------------------------
# Load libraries
# -----------------------------
library(metaboliteIDmapping)
library(AnnotationHub)
library(KEGGREST)
library(dplyr)
library(readxl)
library(writexl)

# -----------------------------
# Step 1: Read metabolite list
# -----------------------------
input_file <- "D:/annotations.xlsx"
met_df <- read_excel(input_file, sheet = "Sheet1")
metabolitesMapping
# -----------------------------
# Step 2: Get HMDB → KEGG mapping
# -----------------------------
mapping_ref <- as.data.frame(metabolitesMapping)

mapped <- merge(
  met_df,
  mapping_ref,
  by.x = "HMDB IDs",
  by.y = "HMDB",
  all.x = TRUE
)


get_compound_mapping <- function(compound_id) {
  # Example compound_id: "C00031" (Glucose)
  tryCatch({
    entry <- keggGet(paste0("cpd:", compound_id))[[1]]
    
    # Extract pathways
    pathways <- if (!is.null(entry$PATHWAY)) {
      data.frame(
        KEGG_Pathway_ID = names(entry$PATHWAY),
        Pathway_Name = entry$PATHWAY,
        stringsAsFactors = FALSE
      )
    } else NULL
    
    # Extract enzymes
    enzymes <- if (!is.null(entry$ENZYME)) {
      data.frame(Enzyme_EC = entry$ENZYME, stringsAsFactors = FALSE)
    } else NULL
    
    # Extract reactions
    reactions <- if (!is.null(entry$REACTION)) {
      data.frame(Reaction_ID = entry$REACTION, stringsAsFactors = FALSE)
    } else NULL
    
    # If no data found, skip
    if (is.null(pathways) && is.null(enzymes) && is.null(reactions)) return(NULL)
    
    # --- Add Pathway Class info ---
    if (!is.null(pathways)) {
      pathway_classes <- sapply(pathways$KEGG_Pathway_ID, function(pid) {
        Sys.sleep(0.5)  # delay to avoid KEGG 403 error
        p_entry <- tryCatch(keggGet(paste0("path:", pid))[[1]], error = function(e) NULL)
        if (!is.null(p_entry$CLASS)) paste(p_entry$CLASS, collapse = "; ") else NA
      })
      pathways$Pathway_Class <- pathway_classes
    }
    
    # Expand all combinations (Compound–Pathway–Enzyme–Reaction)
    df <- expand.grid(
      Compound_ID = compound_id,
      KEGG_Pathway_ID = if (!is.null(pathways)) pathways$KEGG_Pathway_ID else NA,
      Pathway_Name = if (!is.null(pathways)) pathways$Pathway_Name else NA,
      Pathway_Class = if (!is.null(pathways)) pathways$Pathway_Class else NA,
      Enzyme_EC = if (!is.null(enzymes)) enzymes$Enzyme_EC else NA,
      Reaction_ID = if (!is.null(reactions)) reactions$Reaction_ID else NA,
      stringsAsFactors = FALSE
    )
    
    return(df)
    
  }, error = function(e) {
    message(paste("Error retrieving", compound_id, ":", e$message))
    return(NULL)
  })
}
compound_ids <- mapped$KEGG[!is.na(mapped$KEGG)]
compound_ids <- unique(compound_ids)

Sys.sleep_time <- 1.2

compound_mappings <- lapply(compound_ids, function(cid) {
  Sys.sleep(Sys.sleep_time)
  get_compound_mapping(cid)
})

mapping_df <- bind_rows(compound_mappings)

get_enzyme_name <- function(ec) {
  tryCatch({
    entry <- keggGet(paste0("ec:", ec))[[1]]
    entry$NAME[1]
  }, error = function(e) NA)
}

enzyme_names <- unique(mapping_df$Enzyme_EC)
enzyme_map <- data.frame(
  Enzyme_EC = enzyme_names,
  Enzyme_Name = sapply(enzyme_names, get_enzyme_name),
  stringsAsFactors = FALSE
)

mapping_final <- left_join(mapping_df, enzyme_map, by = "Enzyme_EC")

mapping_small <- mapping_final %>%
  select(Compound_ID, Pathway_Name, Enzyme_Name, Pathway_Class) %>%
  distinct()

mapped_small <- mapped %>%
  select('HMDB IDs', 'Original Name', KEGG, Name) %>%
  distinct()

mapping_final2 <- merge(mapped_small, mapping_small, by.x = "KEGG", by.y = "Compound_ID", all.x = TRUE)
write_xlsx(mapping_final2, "D:/metabolite_annotations_simple.xlsx")




# Path to your Excel file
input_file <- "D:/0000.xlsx"

# Read all sheets
disease_metabolite <- read_excel(input_file, sheet = 2)
metabolite_enzyme  <- read_excel(input_file, sheet = 3)
enzyme_pathway     <- read_excel(input_file, sheet = 4)

# Check column names (optional, helps confirm structure)
# names(disease_metabolite)
# names(metabolite_enzyme)
# names(enzyme_pathway)

# Standardize column names to simplify joining
colnames(disease_metabolite) <- c("Disease", "Metabolite")
colnames(metabolite_enzyme)  <- c("Metabolite", "Enzyme")
colnames(enzyme_pathway)     <- c("Enzyme", "Pathway")

# Step 1: Merge disease→metabolite with metabolite→enzyme
disease_enzyme <- merge(disease_metabolite, metabolite_enzyme, by = "Metabolite", all.x = TRUE)

# Step 2: Merge the result with enzyme→pathway
disease_pathway <- merge(disease_enzyme, enzyme_pathway, by = "Enzyme", all.x = TRUE)

# Step 3: Select relevant columns
disease_pathway_final <- disease_pathway %>%
  select(Disease, Pathway) %>%
  distinct() %>%               # remove duplicates
  arrange(Disease, Pathway)

# Step 4: Save output
write_xlsx(disease_pathway_final, "D:/disease_to_pathway.xlsx")


