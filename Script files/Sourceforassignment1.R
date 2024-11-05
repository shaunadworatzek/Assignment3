#Loading needed packages


library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
library("vegan")
library(ggplot2)
library(iNEXT)



#I opened the Nematode data file from bold and explored the variables in the data set 

dfNematodes <- read_tsv(file = "../data/Nematode_data.tsv")

names(dfNematodes)

# I made a new table containing only the variables I needed
dfNematodes.sub <- dfNematodes[, c("processid", "bin_uri","family_name","genus_name", "species_name", "country", "lat", "lon")]

#I made a table with just the families to make sure Onchocercidae was sampled
Nematodes.Families <- table(dfNematodes.sub$family_name)

#Then I made sure there was enough Onchocercidae sampled
names(Nematodes.Families)

hist(x = Nematodes.Families, xlab = "Count of BOLD Records per Family", ylab = "Frequency (No. Families")

sort(table(dfNematodes.sub$family_name), decreasing = TRUE)[1:10]

plot(sort(table(dfNematodes.sub$family_name), decreasing = TRUE)[1:5])

#I made a data frame for just Onchocercidae

dfNematodes.Onchocercidae <- subset(dfNematodes.sub, family_name == "Onchocercidae")

#Then I looked at the number of bins vs species and cleaned up the data by removing missing points

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


#Creating a scatter plot of the Sampling intensity and Species richness at different latitudes

 df.Speciesrichbylat <- df.Onchocercidae %>%
   filter(lat >=0)%>%
   mutate(lat_groups = round(lat)) %>%
   group_by(lat_groups) %>%
   summarize(unique(species_name))%>%
   summarize(total = n())%>%
   as.data.frame()
 
 df.Binbylat <- df.Onchocercidae %>%
   filter(lat >=0)%>%
   mutate(lat_groups = round(lat)) %>%
   group_by(lat_groups) %>%
   summarize(unique(bin_uri))%>%
   summarize(total = n())%>%
   as.data.frame()
 
 df.Samplingbylat <- df.Onchocercidae %>%
   filter(lat >=0)%>%
   mutate(lat_groups = round(lat)) %>%
   group_by(lat_groups) %>%
   summarize(total = n())%>%
   as.data.frame()

 df.graphbylat <- bind_rows(list("Sampling Intensity" = df.Samplingbylat,  "Species richness" = df.Speciesrichbylat, "Bin richness" = df.Binbylat), .id = "id")
 
 gg.graphbylat <- ggplot(df.graphbylat, aes(x= lat_groups, y= total))+
   geom_point(size = 1.5)+
   geom_smooth(method = "lm", color = c("#D55E00"))+
    theme_bw(base_size = 15)+
    ylab("Count")+
    xlab("Latitude")
  
gg.graphbylat + facet_grid(id ~. , scales = "free")

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

#Graph for predicting the total number of species in each latitudinal region with increased sampling 

nem.inext <- iNEXT(listnem, q=0, datatype="abundance", endpoint=150)

ggiNEXT(nem.inext)+
  theme_bw(base_size = 15)+
  scale_color_manual(labels = c("North of 40", "South of 40"), values = c("#D55E00","#0072B2"))+
  scale_shape_manual(labels = c("North of 40", "South of 40"), values = c(19,17))

#Creating a bar graph for the species richness of the species north of 40 degree latitude 

ggplot(dfnemNorth.by.species, aes(x=num_species_north, y=species_name))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Species richness")+
  xlab("Species")



render("Assignment1.Rmd", output_format = "pdf_document")
help(render)
