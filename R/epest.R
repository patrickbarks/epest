
# load packages
library(rgdal)
library(dplyr)
library(tidyr)
library(broom)
library(readr)
library(ggplot2)
library(grid)
library(gganimate)
library(RColorBrewer)

# setwd
setwd('/epest/R')

# read usa counties shapefile
counties <- readOGR('../shp/gz_2010_us_050_00_500k', 'gz_2010_us_050_00_500k')

# remove Alaska (02), Hawaii (15), Puerto Rico (72)
counties <- subset(counties, !STATE %in% c('02', '15', '72'))

# aea projection
namEqualProj4 <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96
                  +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs'

# re-project and remove cols not needed
counties <- counties %>%
  spTransform(CRS(namEqualProj4)) %>% 
  subset(select = c(STATE, COUNTY, NAME, CENSUSAREA))

# create row and fips ids
counties@data$id <- rownames(counties@data)
counties@data$fips <- with(counties@data, paste(STATE, COUNTY, sep = ''))

# tidy and releft_join
countiesTidy <- tidy(counties)
countiesTidy <- left_join(countiesTidy, counties@data, by = 'id')

# read epest dat
epest2009Path <- '../dat/epest_raw/1992-2009/'
epest2012Path <- '../dat/epest_raw/2008-2012/'

epest2009Files <- list.files(epest2009Path)
epest2012Files <- list.files(epest2012Path)

ReadFiles <- function(files, path) read_tsv(paste(path, files, sep = ''))

epest2009List <- lapply(epest2009Files, ReadFiles, path = epest2009Path)
epest2012List <- lapply(epest2012Files, ReadFiles, path = epest2012Path)

epest2009 <- do.call(rbind.data.frame, epest2009List)
epest2012 <- do.call(rbind.data.frame, epest2012List)

# merge state and county fips
epest2009 <- epest2009 %>%
  mutate(STATE_FIPS_CODE = formatC(STATE_FIPS_CODE, width = 2,
                                   flag = '0', mode = 'integer')) %>% 
  unite(fips, STATE_FIPS_CODE, COUNTY_FIPS_CODE, sep = '')

epest2012 <- epest2012 %>%
  mutate(STATE_FIPS_CODE = formatC(STATE_FIPS_CODE, width = 2,
                                   flag = '0', mode = 'integer')) %>% 
  unite(fips, STATE_FIPS_CODE, COUNTY_FIPS_CODE, sep = '')

# any counties not represented in all?
setdiff(epest2009$fips, epest2012$fips)
setdiff(epest2012$fips, epest2009$fips)
setdiff(epest2009$fips, counties@data$fips)
setdiff(epest2012$fips, counties@data$fips)

# correct Dade County (name/fips change in 1997)
# https://www.census.gov/geo/reference/county-changes.html
epest2009$fips[epest2009$fips == 12025] <- 12086

# merge 2009 and 2012 (take 2008-2009 from 2012 file)
intersect(epest2009$YEAR, epest2012$YEAR)
epestMerge <- rbind(filter(epest2009, YEAR < 2008), epest2012)

# merge epest to county areas and scale pest use by county area
epestMergeArea <- left_join(epestMerge, counties@data, by = 'fips') %>%
  mutate(kg_scaled = KG / CENSUSAREA)

# neonicotinoid subset
neonics <- c('CLOTHIANIDIN', 'IMIDACLOPRID', 'ACETAMIPRID', 'THIAMETHOXAM')
epestMergeNeonic <- filter(epestMergeArea, COMPOUND %in% neonics)

# total neonic usage by compound and year
neonicYear <- epestMergeNeonic %>% 
  group_by(COMPOUND, YEAR) %>% 
  summarize(tot_kg = sum(KG))

# plot total neonic usage by year
p1 <- ggplot(neonicYear, aes(x = YEAR, y = tot_kg, col = COMPOUND)) +
  geom_line(size = 2) +
  xlab('Year') + ylab('Total neonicotinoid use (kg)') +
  theme(axis.title = element_text(size = 17),
        text = element_text(size = 15),
        axis.title.x = element_text(margin = margin(0.5, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, 0.5, 0, 0, unit = 'cm')),
        legend.position = c(0.18, 0.8),
        legend.background = element_rect(fill = "grey85"))

#quartz(height = 5, width = 7); print(p1)
ggsave('../img/neonic-yr.tiff', p1, height = 5, width = 7, units = 'in')

# plot scaled neonic usage by year and county
labels <- c('0.00001', '0.001', '0.1', '10')
breaks <- as.numeric(labels)

p2 <- ggplot(epestMergeNeonic, aes(x = YEAR, y = kg_scaled, group = fips)) +
  geom_line(size = 0.2, alpha = 0.05) +
  scale_y_log10(breaks = breaks, labels = labels) +
  xlab('Year') +
  ylab(expression(paste('Neonicotinoid use by county (kg ', mile^-2, ')'))) +
  facet_wrap(~ COMPOUND, ncol = 2) +
  theme(axis.title = element_text(size = 17),
        text = element_text(size = 15),
        axis.title.x = element_text(margin = margin(.5, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .5, 0, 0, unit = 'cm')))

#quartz(height = 5, width = 7); print(p2)
ggsave('../img/neonic-yr-county.tiff', p2, height = 5, width = 7, units = 'in')



## ANIMATION
# My approach was to make a massive tidy df w/ scaled neonic usage by compound,
#  county (incl lat/long), and year, and pass df to gganimate with frame = YEAR.
# Rather than do full join of epest in long format to countiesTidy, I found it
#  much faster to reshape epest to wide format, join to tidy counties shapfile,
#  then reshape back to long format

# make YEAR as.char 
epestMergeNeonic <- epestMergeNeonic %>% 
  mutate(yr = paste('Year', YEAR, sep = '_'))

# reshape to wide format (col for each year)
epestWide <- epestMergeNeonic %>% 
  dplyr::select(fips, COMPOUND, yr, kg_scaled) %>% 
  spread(yr, kg_scaled)

# left_join countiesTidy to epestWide
countiesTidyEpest <- countiesTidy %>% 
  expand(fips, COMPOUND = neonics) %>%
  left_join(countiesTidy, by = 'fips') %>% 
  left_join(epestWide, by = c('fips', 'COMPOUND')) %>%
  gather(yr, kg_scaled, Year_1994:Year_2012)

# get max range of kg_scaled for scale bar
scale_rng <- epestWide %>% 
  dplyr::select(starts_with('Year')) %>%
  as.matrix() %>% 
  range(na.rm = T, finite = T)

# breaks and labels for scale bar
scale_labs <- c('0.00001', '0.001', '0.1', '10')
scale_brks <- as.numeric(scale_labs)

# animate
p3 <- ggplot(countiesTidyEpest,
             aes(x = long, y = lat, group = group, frame = yr)) +
  geom_polygon(aes(fill = kg_scaled)) +
  facet_wrap(~ COMPOUND, ncol = 2) +
  scale_fill_gradient(name = expression(kg / mile^2),
                      na.value = brewer.pal(8, 'Blues')[1],
                      low = brewer.pal(8, 'Blues')[2],
                      high = brewer.pal(8, 'Blues')[8],
                      trans = 'log',
                      limits = scale_rng,
                      breaks = scale_brks,
                      labels = scale_labs) +
  theme(text = element_text(size = 16),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

gg_animate(p3, '../img/neonic-ani.gif', ani.width = 1100, ani.height = 800)
