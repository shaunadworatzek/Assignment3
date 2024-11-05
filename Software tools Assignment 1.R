#Assignment 1
#Alignment; western and eastern monarch butterfly  
# Loading the Biostrings package
library(Biostrings)
# DNA sequences from both west and east monarchs 
seq_west <- DNAString("TATATTTTATTTTTGGAATTTGAGCAGGAATAGTTGGGACATCTTTAAGTCTTTTAATTCGAACAGAATTAGGAACTCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGATTTGGTAATTGATTAGTACCCCTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGACTTTTACCCCCATCATTAATTTTATTAATTTCAAGAAGAATCGTAGAAAATGGTGCAGGAACAGGATGAACAGTTTACCCCCCACTTTCATCAAATATTGCTCATAGAGGATCTTCTGTAGATCTAGCTATTTTTTCTTTACATTTAGCTGGAATTTCATCTATTTTAGGAGCTATTAATTTTATTACTACAATCTTAAATATACGAATTAATAATATAACATTTGATCAAATACCTTTATTTGTTTGAGCAGTAGGTATTACAGCTCTTCTTTTATTACTTTCTTTACCAGTTTTAGCAGGAGCAATTACTATACTTCTTACTGATCGAAATTTAAATACTTCTTTTTTTGATCCTGCTGGTGGAGGAGACCCTATTTTATATCAACATTTATTT")

seq_east <- DNAString("TACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTTGGGACATCTTTAAGTCTTTTAATTCGAACAGAATTAGGAACTCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGATTTGGTAATTGATTAGTACCCCTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGACTTTTACCCCCATCATTAATTTTATTAATTTCAAGAAGAATCGTAGAAAATGGTGCAGGAACAGGATGAACAGTTTACCCCCCACTTTCATCAAATATTGCTCATAGAGGATCTTCTGTAGATCTAGCTATTTTTTCTTTACATTTAGCTGGAATTTCATCTATTTTAGGAGCTATTAATTTTATTACTACAATCTTAAATATACGAATTAATAATATAACATTTGATCAAATACCTTTATTTGTTTGAGCAGTAGGTATTACAGCTCTTCTTTTATTACTTTCTTTACCAGTTTTAGCAGGAGCAATTACTATACTTCTTACTGATCGAAATTTAAATACTTCTTTTTTTGATCCTGCTGGTGGAGGAGACCCTATTTTATATCAACATTTATTT")

# Global alignment; that begins from beginning to the end
global_alignment <- pwalign::pairwiseAlignment(seq_west, seq_east, type = "global")
print(global_alignment)

# Obtaining data from monarch butterflies .tsv file 
# Load necessary libraries
library(stats)
library(tidyverse)
library(viridis)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)

# Load the monarch butterflies data from the .tsv file
monarch_data <- read.delim(file = "./data/denaus.tsv", header = TRUE, sep = "\t")

# Check the structure of the data to understand available columns
str(monarch_data)

# Explore key variables
summary(monarch_data)

# Check for missing values in longitude and latitude columns
sum(is.na(monarch_data$lon))  # Check missing longitude
sum(is.na(monarch_data$lat))  # Check missing latitude

# Remove rows with missing longitude or latitude
monarch_data_clean <- monarch_data[!is.na(monarch_data$lon) & !is.na(monarch_data$lat), ]

# Load natural earth data for country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")


# Ensure the column names ('lon', 'lat') match those in my data
monarch_data_sf <- st_as_sf(monarch_data_clean, coords = c("lon", "lat"), crs = 4326)

# Perform a spatial join to get the country information for each butterfly
joined_data <- st_join(monarch_data_sf, world["name"])


# Visualization 
# Open a PNG device for the whole world plot
png("./data/plot1.png", width = 3000, height = 1800, res = 300)

# Create a world map
ggplot() +
  geom_sf(data = world, fill = "lightblue") +  # Base map
  geom_sf(data = joined_data, aes(color = name), size = 2, alpha = 0.6) +
  labs(title = "Geographical Distribution of Monarch Butterflies Worldwide",
       color = "Country") +
  theme_minimal()
dev.off()

#Exploring North America 
# Load necessary libraries
library(stats)
library(tidyverse)
library(viridis)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)

# Load the monarch butterflies data from the .tsv file
monarch_data <- read.delim(file = "./data/denaus.tsv", header = TRUE, sep = "\t")

# Check the structure of the data to understand available columns
str(monarch_data)

# Explore key variables
summary(monarch_data)

# Display the first few rows of the dataset
head(monarch_data)

# Check for missing values in longitude and latitude columns
sum(is.na(monarch_data$lon))  # Check missing longitude
sum(is.na(monarch_data$lat))  # Check missing latitude

# Remove rows with missing longitude or latitude
monarch_data_clean <- monarch_data[!is.na(monarch_data$lon) & !is.na(monarch_data$lat), ]

# Load natural earth data for country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert monarch data to a spatial object using longitude and latitude
monarch_data_sf <- st_as_sf(monarch_data_clean, coords = c("lon", "lat"), crs = 4326)

# Perform a spatial join to get the country information for each butterfly
joined_data <- st_join(monarch_data_sf, world["name"])

# Visualization for North America
png("./data/monarch_distribution_north_america.png", width = 3000, height = 1800, res = 300)

# Create a map for North America
ggplot() +
  geom_sf(data = world[world$name %in% c("United States of America", "Canada", "Mexico"), ], fill = "lightgray") +  # Base map for North America
  geom_sf(data = joined_data[joined_data$name %in% c("United States of America", "Canada", "Mexico"), ], 
          aes(color = name), size = 1, alpha = 0.6) +
  labs(title = "Geographical Distribution of Monarch Butterflies in North America",
       color = "Country") +
  theme_minimal()
dev.off()


# Filter data for each country
us_data <- joined_data[joined_data$name == "United States of America", ]
canada_data <- joined_data[joined_data$name == "Canada", ]
mexico_data <- joined_data[joined_data$name == "Mexico", ]

# Visualization: Distribution in the United States
png("./data/monarch_distribution_us.png", width = 2500, height = 1800, res = 300)
ggplot() +
  geom_sf(data = world[world$name == "United States of America", ], fill = "lightgray") +  # Base map for the USA
  geom_sf(data = us_data, aes(color = name), size = 3, alpha = 0.6) +
  labs(title = "Geographical Distribution of Monarch Butterflies in the United States of America",
       color = "Country") +
  theme_minimal()
dev.off()

# Visualization: Distribution in Canada
png("./data/monarch_distribution_canada.png", width = 3000, height = 1800, res = 300)
ggplot() +
  geom_sf(data = world[world$name == "Canada", ], fill = "lightblue") +  # Base map for Canada
  geom_sf(data = canada_data, aes(color = name), size = 3, alpha = 0.6) +
  labs(title = "Geographical Distribution of Monarch Butterflies in Canada",
       color = "Country") +
  theme_minimal()
dev.off()

# Visualization: Distribution in Mexico
png("./data/monarch_distribution_mexico.png", width = 3000, height = 1800, res = 300)
ggplot() +
  geom_sf(data = world[world$name == "Mexico", ], fill = "lightblue") +  # Base map for Mexico
  geom_sf(data = mexico_data, aes(color = name), size = 3, alpha = 0.6) +
  labs(title = "Geographical Distribution of Monarch Butterflies in Mexico",
       color = "Country") +
  theme_minimal()
dev.off()


# Add a column to categorize monarch butterflies based on their geographic location (Longitude and Latitude)
monarch_data_clean <- monarch_data_clean %>%
  mutate(region = case_when(
    lon < -100 ~ "Western",   # Longitude less than -100 for Western Monarchs
    lon >= -100 ~ "Eastern"   # Longitude greater than or equal to -100 for Eastern Monarchs
  ))

# Check the first few rows to confirm the categorization
head(monarch_data_clean)

# Visualization: Eastern vs Western Monarchs
png("./data/monarch_east_west_comparison.png", width = 3000, height = 1800, res = 300)
ggplot(monarch_data_clean, aes(x = lon, y = lat, color = region)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~region) +
  labs(
    title = "Comparison of Eastern and Western Monarch Butterflies",
    x = "Longitude",
    y = "Latitude",
    color = "Region"
  ) +
  theme_minimal()
dev.off()

#Novel Question; Why Western monarchs are less than Eastern?
# Obtaining data from monarch butterflies .tsv file 
# Load necessary libraries
library(tidyverse)
# Load the monarch butterflies data from the .tsv file
# Ensure the file path is correct, replace './data/west.tsv' with the actual path if needed
west_data <- read.delim(file = "./data/denaus.tsv", header = TRUE, sep = "\t")


# Filter the data for park and forest habitats, removing rows with missing habitat values
filtered_data <- west_data %>%
  filter(!is.na(habitat) & habitat %in% c("Park", "Forest"))


# Count the number of butterflies in each habitat
habitat_comparison <- filtered_data %>%
  group_by(habitat) %>%
  summarise(total_butterflies = n())  # Count total butterflies for each habitat

# View the comparison
print(habitat_comparison)

# Bar plot to compare total butterfly counts in park and forest
png("./data/monarch_density_north_america.png", width = 3000, height = 1800, res = 300)
ggplot(habitat_comparison, aes(x = habitat, y = total_butterflies, fill = habitat)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Comparison of Western Monarch Butterflies in Parks and Forests",
       x = "Habitat Type",
       y = "Total Butterflies") +
  scale_fill_manual(values = c("blue", "red")) +  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  
dev.off()

