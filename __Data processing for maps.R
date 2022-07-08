# Purpose: UCMR data processing for mapping
# Author: Jahred Liddie

library(tidycensus)
library(tidyverse)
library(sf)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# PFAS data
CWS <- read.csv("../../Merged PFAS data/final_dat_FN_new data.csv")
ucmr <- read.csv("../../Merged PFAS data/UCMR_data.csv")

# ucmr data processing; separate systems serving multiple counties
ucmr.sep <- subset(ucmr, tolower(counties_served) != "not reported") %>% 
  separate_rows(counties_served, sep = ",")

ucmr.sep$state_county <- paste(tolower(ucmr.sep$State), 
                               tolower(ucmr.sep$counties_served), sep = "_")

# text replacements for some municipalities/counties/parishes for easier merging
patterns <- tibble(patterns = c(" municipality", " parish", " municipio",
                                "NA_district of columbia", "ny_new york city",
                                " borough$", " city and$", 
                                " census area$", # alaska county regex
                                "í", "ó", "ü", "ñ", "á"),
                   replacements = c("", "", "", "dc_district of columbia", 
                                    "ny_new york", "", "", "",
                                    "i", "o", "u", "n", "a")
)

# text replacements for UCMR data on municipalities/counties/parishes
# note: U.S. territories are not in ACS, so they will be dropped
for (i in seq_along(patterns$patterns)) {
  
  ucmr.sep$state_county <- str_replace(ucmr.sep$state_county,
                                       pattern = patterns$patterns[i],
                                       replacement = patterns$replacements[i])
  
}

################################################################################
# census data processing for statewide data
demo <- get_acs(geography = "county",
                year = 2018,
                output = "wide",
                variables = c(tot_pop = "B02001_001E",
                              hispanic_tot_pop = "B03001_001E",
                              black_pop = "B02001_003E",
                              hispanic_pop = "B03001_003E",
                              AsAm_pop = "B02001_005E",
                              AmInd_pop = "B02001_004E",
                              white_pop = "B02001_002E",
                              pop_wnh = "B01001H_001E",
                              poverty_pop = "B17001_002E",
                              poverty_tot_pop = "B17001_001E"),
                geometry = TRUE) %>%
  mutate(percBlack = black_pop / tot_pop * 100,
         percHisp = hispanic_pop / hispanic_tot_pop * 100,
         percAsAm = AsAm_pop / tot_pop * 100,
         percAmInd = AmInd_pop / tot_pop * 100,
         percWhite = pop_wnh / tot_pop * 100,
         percColor = (tot_pop - pop_wnh) / tot_pop * 100,
         poverty = poverty_pop / poverty_tot_pop * 100) %>%
  # keeping only estimates, not MoEs
  select(GEOID, starts_with("perc"), tot_pop, poverty,
         NAME)

# 2000-2010 decennial census data
demo_urban <- get_decennial(geography = "county",
                      year = 2010,
                      output = "wide",
                      variables = c(tot_pop = "PCT023001",
                                    tot_pop_urbanq = "H002001",
                                    tot_urban = "H002002",
                                    tot_rural = "H002005"),
                      geometry = FALSE) %>%
  mutate(percurban = tot_urban / tot_pop_urbanq * 100,
         percrural = tot_rural / tot_pop_urbanq * 100) %>%
  # keeping only estimates, not MoEs
  select(GEOID, starts_with("perc"), tot_pop,
         NAME)

# generate FIPS codes and create GEOIDs
fipscodes <- tidycensus::fips_codes

# Create GEOID
fipscodes$GEOID <- paste0(fipscodes$state_code, fipscodes$county_code, sep = "") %>%
  as.integer()

# Create state_county col
fipscodes$state_county <- paste(fipscodes$state, fipscodes$county, sep = '_') %>%
  tolower()

# Clean up for matching
fipscodes$state_county <- str_replace(fipscodes$state_county, " county", "")

# remove everything but GEOID and county name
fipscodes[,c(2, 4:5)] <- NULL

fipscodes$GEOID <- as.character(fipscodes$GEOID)
fipscodes$GEOID <- sprintf("%05s", fipscodes$GEOID)

# merge to get state_county column
demo <- left_join(demo, fipscodes, by = c("GEOID"))
demo_urban <- left_join(demo_urban, fipscodes, by = c("GEOID"))
demo <- left_join(demo, demo_urban %>% select(GEOID, percurban, percrural))

for (i in seq_along(patterns$patterns)) {
  
  demo$state_county <- str_replace(demo$state_county,
                                   pattern = patterns$patterns[i],
                                   replacement = patterns$replacements[i])
  
}

# merge water data with ACS
# first collapse dataset to counties after splitting systems 
# with multiple counties served
CWS_splitcounties <- CWS %>% separate_rows(counties_served, sep = ",")
CWS_splitcounties$state_county <- tolower(paste(CWS_splitcounties$State, 
                                                CWS_splitcounties$counties_served, 
                                                sep = "_"))
CWS_splitcounties <- CWS_splitcounties %>%
  mutate(sources = industries_count + MFTA_count + airport_count)
  
CWS.counties <- CWS_splitcounties %>%
  group_by(state_county) %>%
  summarise(PFOA.avg = mean(PFOA_detect, na.rm = T)*100,
            PFOS.avg = mean(PFOS_detect, na.rm = T)*100,
            PFNA.avg = mean(PFNA_detect, na.rm = T)*100,
            PFBS.avg = mean(PFBS_detect, na.rm = T)*100,
            PFHxS.avg = mean(PFHxS_detect, na.rm = T)*100,
            Any.avg = mean(Any_detect, na.rm = T)*100,
            sources = mean(sources),
            num.CWS = n()
  ) %>%
  ungroup()

ucmr.counties <- ucmr.sep %>%
  group_by(state_county) %>%
  summarise(PFOA.avg = mean(PFOA_detect, na.rm = T)*100,
            PFOS.avg = mean(PFOS_detect, na.rm = T)*100,
            PFNA.avg = mean(PFNA_detect, na.rm = T)*100,
            PFBS.avg = mean(PFBS_detect, na.rm = T)*100,
            PFHxS.avg = mean(PFHxS_detect, na.rm = T)*100,
            PFHpA.avg = mean(PFHpA_detect, na.rm = T)*100,
            Any.avg = mean(Any_detect6, na.rm = T)*100, # detect of 5 compounds
            num.CWS = n()
  ) %>%
  ungroup()

# merge with CWS data
demo.CWS <- left_join(demo, CWS.counties)
demo.CWS <- filter(demo.CWS, !is.na(num.CWS))

demo.ucmr <- left_join(demo, ucmr.counties)
demo.ucmr <- filter(demo.ucmr, !is.na(num.CWS))

st_write(demo.CWS, 
         delete_dsn = TRUE, 
         "../../Merged PFAS data/CWS_formapping.geojson")
st_write(demo.ucmr, 
         delete_dsn = TRUE, 
         "../../Merged PFAS data/UCMR_formapping.geojson")

# also create a ucmr dataset with all CWS merged to census for other analyses
# calculate population-weighted average demographics
ucmr.sep <- left_join(ucmr.sep, demo)

ucmr.sep <- ucmr.sep %>% select(-geometry)

ucmr.CWS <- ucmr.sep %>% 
 group_by(PWSID, PWSName) %>%
 mutate(across(c(starts_with("perc"), "poverty"), ~weighted.mean(.x, tot_pop)),
        tot_pop = sum(tot_pop),
        num_counties = n()) %>%
 ungroup()
 
ucmr.CWS <- ucmr.CWS %>%
  mutate(state_county = ifelse(num_counties > 1, paste(PWSID, "multiple", sep = "_"),
                               state_county))
# name drop!
ucmr.CWS <- ucmr.CWS %>% select(-NAME)

ucmr.CWS <- subset(ucmr.CWS, !duplicated(PWSID, PWSName))

write.csv(ucmr.CWS, "../../Merged PFAS data/ucmr_foranalysis.csv")

