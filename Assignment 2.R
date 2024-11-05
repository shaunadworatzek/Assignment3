# Load necessary libraries
library(Biostrings)      # For DNA sequence manipulation with DNAStringSet
library(ape)             # For phylogenetic analysis and tree visualization
library(phangorn)        # Additional phylogenetic analysis tools
library(dplyr)           # For data manipulation and filtering
library(muscle)          # For sequence alignment with MUSCLE

# Defining file paths for COI and HLA genes
setwd( "../ST assignment 2")
coi_path <- ("../ST assignment 2/COI")
hla_path <- ("../ST assignment 2/HLA")

# List of species for COI and HLA files

species_names <- c("Cat", "Chimpanzee", "GiantPanda", "Gibbon", "Gorilla", 
                   "GrayWolf", "GraySquirrel", "Human", "Mouse", "Rat", 
                   "RedFox", "Rhesus", "Tiger")

# Function to read all FASTA files in a given directory and create a data frame

read_fasta_files <- function(path, species_list, suffix) {
  sequences <- list()
  for (species in species_list) {
    file_path <- file.path(path, paste0(species, suffix, ".fasta"))
    seq_data <- readDNAStringSet(file_path)
    sequences[[species]] <- data.frame(
      processid = species, 
      nucleotides2 = as.character(seq_data),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, sequences)
}

# Read and combine COI and HLA data into data frames

dfCOI <- read_fasta_files(coi_path, species_names, "COI")

dfHLA <- read_fasta_files(hla_path, species_names, "HLA")


# Filter and clean data to retain only the sequences of interest

# Here we make sure each dataframe has unique sequences and processids for consistency

dfCOI <- dfCOI %>% filter(processid %in% species_names)

dfHLA <- dfHLA %>% filter(processid %in% species_names)


# Converting sequences to DNAStringSet format and name them by processid
dfCOI$nucleotides2 <- DNAStringSet(dfCOI$nucleotides2)
names(dfCOI$nucleotides2) <- dfCOI$processid

dfHLA$nucleotides2 <- DNAStringSet(dfHLA$nucleotides2)
names(dfHLA$nucleotides2) <- dfHLA$processid

# Aligning COI sequences using MUSCLE
cat("Aligning COI sequences...\n")
dfCOI.alignment <- DNAStringSet(muscle::muscle(dfCOI$nucleotides2))

# Aligning HLA sequences using MUSCLE
cat("Aligning HLA sequences...\n")
dfHLA.alignment <- DNAStringSet(muscle::muscle(dfHLA$nucleotides2))


# Inspect the alignments
print(dfCOI.alignment)  # Shows aligned COI sequences
print(dfHLA.alignment)  # Shows aligned HLA sequences


# Building Phylogenetic Trees for COI and HLA
# Convert aligned sequences to phyDat format for tree-building

coi_phyDat <- phyDat(as.matrix(dfCOI.alignment), type = "DNA")

hla_phyDat <- phyDat(as.matrix(dfHLA.alignment), type = "DNA")


# Neighbor-Joining Tree for COI

cat("Building Neighbor-Joining Tree for COI...\n")
coi_dist <- dist.ml(coi_phyDat)           # Calculate distance matrix
coi_nj_tree <- NJ(coi_dist)                # Neighbor-Joining tree
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# Maximum Likelihood Tree for COI

cat("Building Maximum Likelihood Tree for COI...\n")
coi_ml_tree <- pml(coi_nj_tree, coi_phyDat)
coi_ml_tree <- optim.pml(coi_ml_tree, model = "HKY")
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Neighbor-Joining Tree for HLA

cat("Building Neighbor-Joining Tree for HLA...\n")
hla_dist <- dist.ml(hla_phyDat)           # Calculate distance matrix
hla_nj_tree <- NJ(hla_dist)                # Neighbor-Joining tree
plot(hla_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

# Maximum Likelihood Tree for HLA

cat("Building Maximum Likelihood Tree for HLA...\n")
hla_ml_tree <- pml(hla_nj_tree, hla_phyDat)
hla_ml_tree <- optim.pml(hla_ml_tree, model = "HKY")
plot(hla_ml_tree$tree, main = "Maximum Likelihood Tree for HLA Sequences")

#Comparing Tree Topologies for COI and HLA
# Here, we calculate a distance metric between the COI and HLA trees

cat("Comparing COI and HLA tree topologies...\n")
tree_distance <- RF.dist(coi_nj_tree, hla_nj_tree)
cat("Robinson-Foulds distance between COI and HLA trees: ", tree_distance, "\n")

# Visualizing Trees Side-by-Side for Comparison
par(mfrow = c(2, 2))  # Set up 2x2 plotting area

# Plot 1: Neighbor-Joining Tree for COI
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# Plot 2: Maximum Likelihood Tree for COI
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Plot 3: Neighbor-Joining Tree for HLA
plot(hla_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

# Plot 4: Maximum Likelihood Tree for HLA
plot(hla_ml_tree$tree, main = "Maximum Likelihood Tree for HLA Sequences")

# Save alignments and tree plots for future reference and analysis
writeXStringSet(dfCOI.alignment, file = "COI_aligned_sequences.fasta")
writeXStringSet(dfHLA.alignment, file = "HLA_aligned_sequences.fasta")
saveRDS(coi_nj_tree, file = "COI_NJ_Tree.rds")
saveRDS(coi_ml_tree, file = "COI_ML_Tree.rds")
saveRDS(hla_nj_tree, file = "HLA_NJ_Tree.rds")
saveRDS(hla_ml_tree, file = "HLA_ML_Tree.rds")

cat("Analysis complete. Results saved.\n")
