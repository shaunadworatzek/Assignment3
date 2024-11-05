# make sure that the working directory is the source file location 


# loading needed packages with some conflicts indicated

library(rentrez)
library(Biostrings)
library(stringr)
library(viridis)
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::rename)
library(rglobi)
library(httr)
library(jsonlite)
library(traitdata)
library(DECIPHER)
library(ape)
library(muscle)
library(phangorn)
library(phytools)
library(ggplot2)
library(ggtree)
library(ggstance)
library(gridExtra)

#### PART 1 - Download data from NCBI cleaning data up for further use----

# investigate the databases on NCBI

entrez_dbs() # Will be using nuccore for my analysis

# investigate the data by searching on the NCBI website and create a web history object to download large amounts of data

Chromadorea18S_search <- entrez_search(db = "nuccore", term = "(Chromadorea[ORGN] AND 18S rRNA AND 200:1500[SLEN])", use_history = T)

Chromadorea18S_search # 13133 sequences were found
Chromadorea18S_search$web_history # A web history object was successfully made

# 18S rRNA was chosen because this is a commonly used marker for identification of nematodes and therefore is variable across but not among species it is also typically used for deep evolutionary relationships

# loading source function to download larger data sets from NCBI

source("Entrez_Functions.R")

# fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Chromadorea[ORGN] AND 18S rRNA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Chromadorea_18S")

# combing the many 100 Chromadorea sequence files into one data frame and moving it into R

Chromadorea_18S <- MergeFastaFiles(filePattern = "Chromadorea_18S*")

# cleaning up sequence data to only have the sequence and the species name by separating the first column into different groups and removing any missing sequence or species data and any species with mismatched data

df_Chromadorea_18S <- Chromadorea_18S %>%
  separate(Title,
    into = c("Code", "Genus", "Species", "rest"),
    sep = " ",
    extra = "merge",
    fill = "right",
    convert = FALSE
  ) %>%
  mutate(Species = paste(Genus, Species)) %>%
  mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>%
  select(Code, Species, Sequence, spaces) %>%
  filter(spaces == 1, !is.na(spaces), !is.na(Species)) %>%
  dplyr::rename(Parasite_species = Species)

length(unique(df_Chromadorea_18S$Parasite_species)) # 1289 Chromadorea species are remaining in the data

# creating a new column for sequence length to make sure that the sequences are all of a similar length
df_Chromadorea_18S <- df_Chromadorea_18S %>%
  mutate(seqlength = nchar(df_Chromadorea_18S$Sequence))

# visualizing the sequence lengths with a histogram

hist(x = df_Chromadorea_18S$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")

# selecting for sequence lengths between 500 and 1000 base pairs to reduce the number of sequences in analysis and therefore reduce computational time

df_Chromadorea_18S <- df_Chromadorea_18S %>%
  filter(seqlength <= 1000) %>%
  filter(seqlength >= 500)

# visualizing the new sequence lengths with a histogram

hist(x = df_Chromadorea_18S$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")

length(unique(df_Chromadorea_18S$Parasite_species)) # 1042 species are present after filtering out the larger sequences, which is okay for this project because there were too many species anyway

#### PART 2 -Downloading trait data on parasitism from GloBI using pagination----


# using pragnation to download all the available data instead of a smaller subset of data

# defining the limit and skip for pagination

limit <- 100
skip <- 0

# defining the source taxon and interaction type
source_taxon <- "Chromadorea"
interaction_type <- "parasiteOf"

# creating an empty data frame to store all of the results
Chromadorea_trait <- data.frame()

# making a loop to fetch all data through from the global interactions database
repeat {
  # building the URL with query parameters
  url <- paste0(
    "https://api.globalbioticinteractions.org/interaction?",
    "sourceTaxon=", source_taxon,
    "&interactionType=", interaction_type,
    "&limit=", limit,
    "&skip=", skip
  )

  # send a GET request to the API
  response <- GET(url)

  # parse the response as JSON
  interactions_chunk <- fromJSON(content(response, "text"), flatten = TRUE)

  # if no results are retrieved, break the loop
  if (length(interactions_chunk$data) == 0) {
    break
  }

  # add the results to the created data frame
  Chromadorea_trait <- rbind(Chromadorea_trait, as.data.frame(interactions_chunk$data))

  # update the skip value for pagination
  skip <- skip + limit
}

# looking at the column names

names(Chromadorea_trait) # they are not defined and are currently V1...

# investigating the data and removing any missing data or species of both the parasite and host which are not in the correct format, then selecting for only the distinct taxa so that each species has its own entery

length(unique(Chromadorea_trait$V2)) # 4215 Chromadorea species are present before cleaning up the data set

df.Chromadorea_trait <- Chromadorea_trait %>%
  mutate(spaces_source = str_count(V2, "[\\s\\.\\d]"), spaces_target = str_count(V12, "[\\s\\.\\d]")) %>%
  filter(spaces_source == 1, !is.na(spaces_source), spaces_target == 1, !is.na(spaces_target)) %>%
  select(Species = V2, Interaction = V10, Target_species = V12, Target_taxonomy = V13) %>%
  distinct()

length(unique(df.Chromadorea_trait$Species)) # 3395 Chromadorea species are present after cleaning up the data set


#### PART 3 - Downloading data from the ecological trait database on the data set PANtheria for mammalian orders and adult body mass data  ----

# opening Mammal trait data into R and selecting for order, body size, genus and species, then combing genus and species into one column called target species

mammal.trait <- pantheria %>%
  select(AdultBodyMass_g, Order, Genus, Species) %>%
  filter(!is.na(AdultBodyMass_g)) %>%
  unite(Target_species, Genus, Species, sep = " ")

# conformation of a suitable number of orders for analysis as in not too many and not too little (if too little then I would move down to family)

length(unique(mammal.trait$Order)) # 29 Orders were present so this is suitable moving forward

# combing mammal trait data with global interactions data

df.traitdata <- full_join(mammal.trait, df.Chromadorea_trait, join_by("Target_species" == "Target_species"))

# looking at the biomass data

summary(df.traitdata$AdultBodyMass_g) # The min is 2 and the max is 154321304, there are 23618 NAs

# cleaning up the combined data so that all the missing data is removed, only the Parasite species, host order, and host body size are present, filtering out species which have exceptionally large body mass (its a whale) as although this is interesting the magnitude difference makes it difficult to analyze the other data

df.Chromadorea_trait <- df.traitdata %>%
  filter(!is.na(Species), !is.na(AdultBodyMass_g)) %>%
  select(Species, AdultBodyMass_g, Order) %>%
  distinct(Species, .keep_all = TRUE) %>%
  filter(AdultBodyMass_g <= 1000000) %>%
  dplyr::rename(Parasite_species = Species, Mammillian_order = Order)

# checking on the state of the combined data

view(df.Chromadorea_trait)


length(unique(df.Chromadorea_trait$Mammillian_order)) # 22 orders in data only 7 were lost when combing with parasites
length(unique(df.Chromadorea_trait$Parasite_species)) # 1269 Nematode species are present in the data set


#### PART 4 - preparing the data for determing the centroid sequence for each species----

# combing trait data with the sequence data to reduce the number of species in the centroid calculation and therefore reduce the computation time (since all the species in the two data sets are not the same I might as well mimit the number of species being used in further analysis now)

df.Chromadorea_all <- full_join(df.Chromadorea_trait, df_Chromadorea_18S, join_by("Parasite_species" == "Parasite_species"))

# removing any missing data

df.Chromadorea_all <- df.Chromadorea_all %>%
  filter(!is.na(Sequence), !is.na(AdultBodyMass_g))

# looking at the combined data and the number of species of parasites remaining

view(df.Chromadorea_all)
length(unique(df.Chromadorea_all$Parasite_species)) # 141 unique Chromadorea species when the sequence data and the trait data are combined
unique(df.Chromadorea_all$Parasite_species)

# Making the sequences into a DNAstringset

class(df.Chromadorea_all$Sequence)
df.Chromadorea_all$Sequence2 <- DNAStringSet(df.Chromadorea_all$Sequence)
class(df.Chromadorea_all$Sequence2)

# Looking at the sequences on an online browser

BrowseSeqs(df.Chromadorea_all$Sequence2)

#### - PART 5 - determing the centroid sequence ----

# Grouping sequences by species

grouped_sequences <- split(df.Chromadorea_all$Sequence2, df.Chromadorea_all$Parasite_species)

# Creating the function to find the centroid

calculate_centroid <- function(seqs) {
  # performing multiple sequence alignment
  alignment <- DNAStringSet(muscle::muscle(seqs, gapOpening = -3000))

  # converting the alignment to a DNAbin object to calculate distances
  alignment_dnabin <- as.DNAbin(alignment)

  # calculating pairwise distance matrix using TN93 model
  dist_matrix <- dist.dna(as.DNAbin(alignment), model = "TN93")

  # calculating centroid (sequence with the lowest sum of pairwise distances)
  centroid_index <- which.min(rowSums(as.matrix(dist_matrix)))

  # extracting the sequence of the centroid
  centroid_sequence <- as.character(alignment[centroid_index])

  return(centroid_sequence)
}

# Apply the centroid calculation for each specimens and save it as a data frame

centroids <- lapply(grouped_sequences, calculate_centroid)
class(centroids)

#### - PART 6 - phylogeny tree construction with parsimony and maximum likely and comparison between the two ----

# separating each nucleotide into its own column and keeping the species names as the row names and then converting the data frame into a matrix

df.centroids <- centroids %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(Sequence_split = strsplit(as.character(V1), "")) %>%
  unnest_wider(Sequence_split, names_sep = "V1_") %>%
  mutate(rowname = gsub("\\.", " ", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  select(-V1) %>%
  as.matrix()

class(df.centroids) # is a matrix/array

# matrix centroid conversion to a phyDat format

DNA.Chromadorea_centroid <- as.phyDat(df.centroids)

class(DNA.Chromadorea_centroid) # is a "phyDat" object
length(DNA.Chromadorea_centroid) # has the 141 species as seen before

# creating a distance matrix

dist.Chromadorea <- dist.ml(DNA.Chromadorea_centroid, ratio = TRUE, model = "F81") # no this is definatly the one I should use

# creating a tree using the neighbor joining method

tree.Chromadorea <- NJ(dist.Chromadorea)

# determining the parsimony score of the tree

parsimony(tree.Chromadorea, DNA.Chromadorea_centroid, method = "fitch")
parsimony(tree.Chromadorea, DNA.Chromadorea_centroid, method = "sankoff") # 28363 is the parsimony score using either the fitch or sankoff method

# finding the best tree (or list of trees)

treeRatchet <- pratchet(DNA.Chromadorea_centroid, start = tree.Chromadorea, maxit = 100, minit = 5, k = 5, trace = 0)

# returns a tree with edge lengths according to the ACCTRAN criterion

tree.best <- acctran(treeRatchet, DNA.Chromadorea_centroid)

# checking the bootstrapping values

plotBS(midpoint(tree.best), type = "phylogram", cex = 0.6, show.tip.label = FALSE)

# model testing to find the best model

modelTest <- modelTest(DNA.Chromadorea_centroid, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # F81+G(4) is the best model because it has the lowest AIC, second lowest BIC and high loglik
summary(modelTest)
view(modelTest)

# making a tree with maximum likelihood using the F81+G(4) method

fit_mt <- pml_bb(DNA.Chromadorea_centroid, model = "F81+G(4)", control = pml.control(trace = 0))

# finding the bootstrap values for the maximum likelihood tree

bs <- bootstrap.pml(fit_mt, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

# checking the bootstrpaing values

plotBS(midpoint(fit_mt$tree), type = "phylogram", cex = 0.6, show.tip.label = FALSE) # bootstrap values seem to be

# comparing the bootstrap values

# bootstrap values for tree.best
print(tree.best$node.label)

average.tree.best <- mean(as.numeric(tree.best$node.label), na.rm = TRUE)

print(average.tree.best) # 0.4883 is the average bootstrap value

# boostrap values for fit_mt$tree
print(fit_mt$tree$node.label)

average.tree.mt <- mean(as.numeric(fit_mt$tree$node.label), na.rm = TRUE)

print(average.tree.mt) # 0.5709 is the average bootstrapping value for the maximum likelihood tree

# meaning that there is more support for the maximum liklihood tree but this support still is only moderate

# comparing different trees

RF.dist(tree.best, fit_mt$tree) # 196 and the maximum value for unrooted trees is (2×(140−1)=279) which indicates a relativley large discrepancy between trees

# Making a plot comparing the two different trees

# creating a ggtree for the parsimonous tree

parsimony <- ggtree(tree.best) + ggtitle("Parsimony")

# creating a ggtree for the maximum liklihood tree

maxliklihood <- ggtree(fit_mt$tree) + ggtitle("Maximum liklihood")

# combining the two plots

grid.arrange(parsimony, maxliklihood, ncol = 2)

# Figure 1. Comparison of phylogenies of parasitic nematodes with a mammalian host computed either based on the most parsimonious tree or the maximum likelihood.

#### - PART 7 - ggtree construction for the orders of the mammalian host represented at the tips  ----

# making the centroids into a data frame so that the species can be ordered in the same manor as they are in the trait data

df.DNA <- centroids %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(rowname = gsub("\\.", " ", rowname)) %>%
  dplyr::rename(Parasite_species = rowname)

# combing the centroid data frame from above with the trait data from part 3 to make sure the trait data only contains the species in the tree

df.Trait_DNA <- left_join(df.DNA, df.Chromadorea_trait, join_by("Parasite_species" == "Parasite_species"))

# filtering the trait and centroid data to select for the parasite species, mammalian order, Adult body mass, and then creating a new column for the body mass groups

extra.info <- df.Trait_DNA %>%
  select(Parasite_species, Mammillian_order, AdultBodyMass_g) %>%
  mutate(Mass_grouping = cut(df.Trait_DNA$AdultBodyMass_g, breaks = 5))


length(extra.info$Parasite_species) # 141 species confirming the length is correct
length(unique(extra.info$Mammillian_order)) # 11 mammalian orders are present
class(extra.info) # data frame

# the next few steps are done to make sure that the trait data and the tree have the same order of parasitic species for the ggtree construction

# getting the tip labels from the tree
tree_tips <- fit_mt$tree$tip.label

# getting the species names from the trait data
trait_species <- extra.info$Parasite_species

# checking if all species names are present in the tree
all(trait_species %in% tree_tips)

# reordering the trait data to match the order of the tree tip labels
extra.info <- extra.info[match(tree_tips, extra.info$Parasite_species), ]

node <- 1:Ntip(tree.best)

# making the tree in a circular layout with the tips having points representing the mammilian order of the host of each particular species

ggtree(fit_mt$tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(color = extra.info$Mammillian_order[node]), size = 2) +
  scale_color_manual(values = c("#76EEC6", "#CD3333", "#7AC5CD", "#D2691E", "#EE6A50", "#2F4F4F", "#698B69", "#EEC900", "#EE6AA7", "#C1CDC1", "#CD5C5C"), "Mammalian Order") +
  theme_tree2(legend.position = "right")

# Figure 2. Parasitic nematode phylogeny with the order of the mammalian host represented on the tree tips.

#### PART 8 - Testing for phylogenetic conservatism and body size plus ggtree with bargrpah of body size ----

# creating a vector for the adult body mass

traitvector_Chromadoera <- setNames(extra.info$AdultBodyMass_g, extra.info$Parasite_species)

# making sure the species order is the same in both

traitvector_Chromadoera <- traitvector_Chromadoera[tree.best$tip.label]

# estimating the lambda parameter for host body size with the maximum likelihood method

lambda.mt <- phylosig(fit_mt$tree, traitvector_Chromadoera, method = "lambda", test = TRUE, nsim = 1000, se = NULL, start = NULL)

# lambda is equal to 7.331374e-05 which means that the traits are not conserved

# estimating the lambda parameter for host body size with the parsimony

lambda.ps <- phylosig(tree.best, traitvector_Chromadoera, method = "lambda", test = TRUE, nsim = 1000, se = NULL, start = NULL)

# lambda is equal to 0.1284 which shows some slight conservation of traits

# creating a ggtree object from the tree in the form of a cladogram (this is how it differs from before) for both max likelihood and most parsimonious

cladogram.pars <- ggtree(tree.best, branch.length = "none")

cladogram.max <- ggtree(fit_mt$tree, branch.length = "none")

# using faceting to combine the tree with a bar graph for body size

facet_plot(cladogram.pars, panel = "Adult body mass (g)", geom = geom_barh, data = extra.info, mapping = aes(x = AdultBodyMass_g, fill = Mass_grouping), stat = "identity") +
  scale_fill_manual(
    name = "Mass grouping", labels = c("2 to 135000", "135000 to 270000", "270000 to 406000", "406000 to 541000", "541000 to 677000"),
    values = c("#2F4F4F", "#CD3333", "#7AC5CD", "#D2691E", "#EE6A50")
  )

facet_plot(cladogram.max, panel = "Adult body mass (g)", geom = geom_barh, data = extra.info, mapping = aes(x = AdultBodyMass_g, fill = Mass_grouping), stat = "identity") +
  scale_fill_manual(
    name = "Mass grouping", labels = c("2 to 135000", "135000 to 270000", "270000 to 406000", "406000 to 541000", "541000 to 677000"),
    values = c("#2F4F4F", "#CD3333", "#7AC5CD", "#D2691E", "#EE6A50")
  )

# Figure 3.Plot of the phylogeny of parasitic nematodes with the body mass of their respective host shown on the left the colours represent the different groups of weight in grams where A) is the most parsimonious and B) is thr max liklihood
