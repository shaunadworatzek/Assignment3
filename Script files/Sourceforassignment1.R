## Loading required packages for data analysis and processing.

library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
library("vegan")
library(ggplot2)
library(iNEXT)
##New libraries 
library(furrr) 
plan(multisession)  # Set up parallel processing
library(maps)
library(sf)
library(plotly)

## Opened the Nematode dataset from BOLD and conducted an initial exploration of the variables

dfNematodes <- read_tsv(file = "./data/Nematode_data.tsv")

names(dfNematodes)

## Created a new table containing only the essential variables required for analysis
dfNematodes.sub <- dfNematodes[, c("processid", "bin_uri","family_name","genus_name", "species_name", "country", "lat", "lon")]

## Created a table with only the family-level data to confirm the presence of Onchocercidae in the sample
Nematodes.Families <- table(dfNematodes.sub$family_name)

## Generated a data frame filtered specifically for the Onchocercidae family
dfNematodes.Onchocercidae <- subset(dfNematodes.sub, family_name == "Onchocercidae")


## The following code  will ensure that you are working with a clean data set without missing values.

# Function to log rows with missing data
log_missing_data <- function(data, file = "missing_data_log.csv") {
  missing_data <- data[!complete.cases(data), ]
  if (nrow(missing_data) > 0) {
    write.csv(missing_data, file, row.names = FALSE)
    cat("Logged", nrow(missing_data), "rows with missing values to", file, "\n")
  }
}

# Check for missing values in key columns and filter them out
required_columns <- c("species_name", "bin_uri", "lat", "country")

df.Onchocercidae <- dfNematodes.Onchocercidae %>%
  filter(complete.cases(select(., all_of(required_columns))))

# Inform if any records were removed due to missing data
initial_count <- nrow(dfNematodes.Onchocercidae)
filtered_count <- nrow(df.Onchocercidae)

if (initial_count > filtered_count) {
  cat("Removed", initial_count - filtered_count, "rows with missing values in key columns.\n")
}

# Log missing values for troubleshooting
log_missing_data(dfNematodes.Onchocercidae)

## Verified that the sample size for the Onchocercidae family was sufficient for analysis
names(Nematodes.Families)

hist(x = Nematodes.Families, xlab = "Count of BOLD Records per Family", ylab = "Frequency (No. Families")

sort(table(dfNematodes.sub$family_name), decreasing = TRUE)[1:10]

plot(sort(table(dfNematodes.sub$family_name), decreasing = TRUE)[1:5])


## Examined the relationship between the number of bins and species, then cleaned the data by removing any missing values.

length(unique(dfNematodes.Onchocercidae$bin_uri))
length(unique(dfNematodes.Onchocercidae$species_name))

dfNematodes.Onchocercidae$spaces <- str_count(string = dfNematodes.Onchocercidae$species_name, pattern = "[\\s\\.\\d]")

df.Onchocercidae<- subset(dfNematodes.Onchocercidae, spaces == 1 & is.na(spaces) == F & is.na(species_name) == F & is.na(bin_uri) == F)
length(unique(df.Onchocercidae$species_name))
length(unique(dfNematodes.Onchocercidae$bin_uri))
rm(df.Onchocercidae)

df.Onchocercidae<- subset(dfNematodes.Onchocercidae, spaces == 1 & is.na(spaces) == F & is.na(species_name) == F & is.na(bin_uri) == F & is.na(country) == F & is.na(lat) == F)
length(unique(df.Onchocercidae$species_name))
length(unique(dfNematodes.Onchocercidae$bin_uri))


# Creating a scatter plot of the Sampling intensity and Species richness at different latitudes
# Improving summarize_by_latitude Function with Error Handling 
# Function to summarize data by latitude groups
summarize_by_latitude <- function(data, column_name) {
  if (!column_name %in% colnames(data)) {
    stop("Column not found in the data frame.")
  }
  data %>%
    filter(lat >= 0) %>%
    mutate(lat_groups = round(lat)) %>%
    group_by(lat_groups) %>%
    summarize(total = n_distinct(.data[[column_name]])) %>%
    as.data.frame()
}

## Parallelize summarization by latitude
results <- future_map(
  list("species_name", "bin_uri", "lat"),
  ~ summarize_by_latitude(df.Onchocercidae, .x)
)

# Extract results for each summary
df.Speciesrichbylat <- results[[1]]
df.Binbylat <- results[[2]]
df.Samplingbylat <- results[[3]]

# Combine data frames for plotting
df.graphbylat <- bind_rows(
  list(
    "Sampling Intensity" = df.Samplingbylat,
    "Species richness" = df.Speciesrichbylat,
    "Bin richness" = df.Binbylat
  ),
  .id = "id"
)

## Enhanced scatter plot with improved color scheme and labels

gg.graphbylat <- ggplot(df.graphbylat, aes(x= lat_groups, y= total)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", color = "#0072B2") +
  theme_bw(base_size = 15) +
  labs(
    title = "Sampling Intensity and Species Richness Across Latitudes",
    x = "Latitude (Rounded)",
    y = "Count"
  ) +
  facet_grid(id ~ ., scales = "free")

print(gg.graphbylat)


#Exploring by latitude data sets and finding the statistical results of the regression 

hist(df.Binbylat$total)
hist(df.Speciesrichbylat$total)
hist(df.Samplingbylat$total)


modelbinrichness <- lm(total ~ lat_groups, data = df.Binbylat)
modelspeciesrichness <- lm(total ~ lat_groups, data = df.Speciesrichbylat)
modelsampling <- lm(total ~ lat_groups, data = df.Samplingbylat)

summary(modelbinrichness)
summary(modelspeciesrichness)
summary(modelsampling)

# Statistical tests for latitude groups using ANOVA
anova_binrichness <- aov(total ~ lat_groups, data = df.Binbylat)
anova_speciesrichness <- aov(total ~ lat_groups, data = df.Speciesrichbylat)
anova_sampling <- aov(total ~ lat_groups, data = df.Samplingbylat)

# Display the ANOVA summary for each test
cat("ANOVA Results for Bin Richness by Latitude:\n")
print(summary(anova_binrichness))

cat("\nANOVA Results for Species Richness by Latitude:\n")
print(summary(anova_speciesrichness))

cat("\nANOVA Results for Sampling Intensity by Latitude:\n")
print(summary(anova_sampling))

#Making two species accumulation curve, one for above 40 lat, and one between 0 and 40 lat with iNEXT to incorporate extrapolation 
#Creating a new data set for only between 40 degrees latitude and the equator

df.nemsouth <- df.Onchocercidae %>%
  filter(lat <=40, lat >=0)

#creating a new data set for only above 40 degrees latitude

df.nemNorth <- df.Onchocercidae %>%
  filter(lat >= 40)

#Grouping per species and counting the number of each species in both the north and south subset 

dfnemsouth.by.species <- df.nemsouth %>%
  group_by(species_name) %>%
  count(species_name)

colnames(dfnemsouth.by.species)[2] <- "num_species_south"

dfnemNorth.by.species <- df.nemNorth %>%
  group_by(species_name) %>%
  count(species_name)

colnames(dfnemNorth.by.species)[2] <- "num_species_north"

#Converting the number of species data into one list

num_species_northd <- as.double(as.character(dfnemNorth.by.species$num_species_north))

num_species_southd <- as.double(as.character(dfnemsouth.by.species$num_species_south))

listnem <- list(num_species_northd, num_species_southd)

## Graph for predicting the total number of species in each latitudinal region with increased sampling 
## I made changes in order to improvement
nem.inext <- iNEXT(listnem, q=0, datatype="abundance", endpoint=150)

ggiNEXT(nem.inext) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Species Accumulation Curve by Latitudinal Region",
    x = "Sample Size",
    y = "Species Richness",
    color = "Region",
    shape = "Region"
  ) +
  scale_color_manual(
    values = c("#0072B2", "#D55E00"),
    labels = c("North of 40°", "South of 40°")
  ) +
  scale_shape_manual(
    values = c(19, 17),
    labels = c("North of 40°", "South of 40°")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "top", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5),
    shape = guide_legend(title.position = "top", title.hjust = 0.5)
  )

# Creating a bar graph for the species richness of the species north of 40 degree latitude 

ggplot(dfnemNorth.by.species, aes(x=num_species_north, y=species_name))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Species richness")+
  xlab("Species")

# Geospatial visualization of sampling locations
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = dfNematodes.sub, aes(x = lon, y = lat), color = "blue", size = 1) +
  labs(title = "Sampling Locations for Nematodes", x = "Longitude", y = "Latitude") +
  theme_minimal()



