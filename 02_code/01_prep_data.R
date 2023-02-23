## 01_prep_data.R
## Run pubmed queries with package, clean and standardize metadata for analysis
rm(list = ls())

## Load packages
library(pacman)
pacman::p_load(tidyverse, rio, data.table, plyr, janitor, assertable, lubridate, furniture, zoo, viridis, classInt)

## Import
files <- data.table(paths = list.files("01_data", pattern = "txt"))
files[, c("organism", "datestamp") := tstrsplit(paths, split = "_")][, datestamp := as.integer(str_remove(datestamp, "\\.txt"))]
files <- files[, .SD[datestamp == max(datestamp)], by = .(organism)]$paths

df <- import_files(files, folder = "01_data", FUN = function(f){
  
  d <- fread(f) %>% clean_names()
  data_organism <- str_remove(f, "01_data\\/") %>% tstrsplit(., split = "_") %>% .[[1]]
  d[, organism := data_organism]
  
  return(d)
  
})

## Pathogen metadata
search_df <- import("01_data/fpp_search_builder.xlsx") %>% as.data.table

## Filter to only the variables we care about: drop, AD a little too messy about including country codes, PT, keyword vars for now to keep prep simple MH OT
df[, v2 := NULL]
df <- df[v1 %in% c("PMID", "DP", "TI", "AB", "AU", "TA", "PL", "LA")]

## Indexing by PMID to make grouping easier
df[v1 == "PMID", pmid_count := seq(.N), by = v1]
df[, pmid_count := na.locf(pmid_count)]
pmid_map <- unique(df[v1 == "PMID", .(PMID = v3, pmid_count)])
df <- merge(df, pmid_map, by = "pmid_count", all.x = T)
df <- df[v1 != "PMID"]
df[, pmid_count := NULL]

## THis is weird just fix it
df[, v4:= str_replace_na(v4, "")][, v3 := str_c(v3,v4)][, v4 := NULL]

## Pulling out first author and associated information, dropping rest
df[v1 %in% c("AU", "AD"), index := seq(.N), by = .(v1, PMID, organism)]
df <- df[is.na(index)|index == 1]
df[, index := NULL]

## Where languages are multiple, if english exists keep only english
#df[v1 == "LA", n_lang := .N, by = .(PMID, organism)]
df[, lang_eng := any(v3 == "eng"), by = PMID]
df[lang_eng == T & v1 == "LA", v3 := "eng"]
df[, lang_eng := NULL]
df <- unique(df)

## Casting out wide by citation
df <- dcast(df, PMID + organism ~ v1, value.var = "v3")


## turning pubdate into year
df[, DP := as.integer(str_extract(DP, "^[:digit:]{4}"))]


## Extract country from first author affiliation: later

## Clean up column names
setnames(df, c("DP", "TI", "AB", "AU", "TA", "PL", "LA"), c("pub_year", "title", "abstract", "first_author", "journal", "pub_location", "pub_lang"))
df <- clean_names(df)

## Binning years
df[, pub_year_window_5 := round_any(pub_year, 5, floor)]

## Merging on fpp groups
df <- merge(df, search_df[, .(organism = fpp, priority)], by = "organism", all.x = T)

## Unit of analysis is organism-pmid (only one where this seems sketchy is scedosporium spp vs l prolificans)

## Exploratory data analysis


## Wishlist


## Tables
## Priority -> Organism
## Number of citations, citations in 2022, unique journals, most frequent journal, unique languages, most frequent non-english lang
n_df <- df[, .N, by = .(organism, priority)]
n22_df <- df[pub_year == 2022, .N, by = organism] %>% setnames(., "N", "N (2022)")
njournal_df <- unique(df[, .(organism, journal)])[, .N, by = organism] %>% setnames(., "N", "Unique Journals")
nlang_df <- unique(df[, .(organism, pub_lang)])[, .N, by = organism] %>% setnames(., "N", "Unique Languages")

mlang_df <- df[, .N, by = .(organism, pub_lang)][order(-N)][, .SD[N != max(N)], by = organism][, .SD[N == max(N)], by = organism]
mlang_df <- mlang_df[, .(pub_lang = paste(pub_lang, collapse = ", ")), by = .(organism)]  %>% setnames(., "pub_lang", "Top non-English")

mjournal_df <- df[, .N, by = .(organism, journal)][order(-N)]
mjournal_df <- mjournal_df[, .(journal = paste(journal, collapse = ", ")), by = .(organism, N)]
mjournal_df[, rank := seq(.N), by = organism]
mjournal_df <- mjournal_df[rank <= 3]
mjournal_df[, journal := sprintf("(%d) %s", rank, journal), by = organism]
mjournal_df <- mjournal_df[, .(journal = paste(journal, collapse = ", ")), by = .(organism)]  %>% setnames(., "journal", "Top Journals")

t1 <- list(n_df, n22_df, njournal_df, mjournal_df, nlang_df, mlang_df) %>% reduce(., function(x, y) merge(x, y, by = "organism"))
t1 <- t1[order(-priority, -organism)]
write.csv(t1, "03_outputs/t1.csv", row.names = F)

## Figs
## Data density heatmap over time
dd_df <- df[, .N, by = .(organism, priority, pub_year_window_5)]
jbreaks <- classIntervals(dd_df$N, style = "quantile", n = 10)
dd_df[, N_discrete := cut(N, breaks = c(jbreaks$brks), include.lowest = T)]
dd_df <- dd_df[order(-priority, -organism)]
dd_df[, organism := factor(organism, levels = unique(dd_df$organism))]
hm <- ggplot(data = dd_df[pub_year_window_5 >= 1990], aes(x = pub_year_window_5, y = organism)) + 
  geom_tile(aes(fill = N_discrete), color = "black") + 
  scale_fill_brewer(palette = "RdYlGn", labels = c("Fewer", rep("", 8), "More")) + theme_minimal() + labs(fill = "Unique citations")

pdf("03_outputs/hm.pdf", height = 6, width = 4.5)
print(hm)
dev.off()

## Top journal heatmap

## Lolipop/Bar: country/region
## Map (less critical

pl_df <- df[, .N, by = .(pub_location, organism, priority)][order(-N)]
pl_df[, pub_location := factor(pub_location, levels = unique(pl_df$pub_location))]
pl_df[, organism := factor(organism, levels = unique(pl_df$organism))]
pl_bar <- ggplot(data = pl_df[N > 20], aes(x = pub_location, y = N)) + geom_bar(stat = "identity", aes(fill = priority)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")

pdf("03_outputs/pl_bar.pdf", height = 6, width = 4.5)
print(pl_bar)
dev.off()

## Histogram/heatmap over time
pr_hist_df <- df[, .N, by = .(priority, pub_year_window_5)]
hist <- ggplot(data = pr_hist_df, aes(x = pub_year_window_5, y = N, color = priority, fill = priority)) + geom_bar(stat = "identity", position = "stack") + theme_minimal()


## Map/Barchart of pub locations: by WHO region

## Frequency by journal (top 10 - all fungal publications, look at impact factor, top journal by organism)
journ_df <- df[, .N, by = .(journal, priority)][order(priority, -N)]
journ_org_df <- df[, .N, by = .(journal, organism)][order(organism, -N)]

# most common journal by priority tier, organism
journ_df[, .SD[N == max(N)], by = priority]
journ_exp <- journ_org_df[, .SD[N == max(N)], by = organism]

## top 10 journals by priority tier
journ_df[, .SD[1:10,], by = .(priority)][, unique(journal)]

## Country map by organism


## Keeping only study types that we want. Can do this later
## Determining whether study type is review or meta analysis
## Dropping study types we want to drop: retracted, letters, newspaper, clinical trial protcol

## Care about
## Glossary here: https://www.nlm.nih.gov/bsd/mms/medlineelements.html
## PL = place of publication
## PMID, DP, TI, AB, AU/AD (first instance = first author), PT keep all, PL (publicaiton location?), TA, MH (mesh), OT (other terms)

## Can maybe see if we can merge impact factor on to journal and compare with other conditions



## Do shit
## Sample query
q <- '(aspergill*[Title/Abstract] OR neuroaspergill*[Title/Abstract] OR fumigatus[Title/Abstract]) AND (incidence[Title/Abstract] OR prevalence[Title/Abstract] OR mortality[Title/Abstract] OR case fatality[Title/Abstract] OR identification[Title/Abstract] OR outbreak[Title/Abstract] OR cluster[Title/Abstract] OR burden[Title/Abstract]) AND 
("Epidemiologic Studies"[Mesh] OR "Population Surveillance"[Mesh] OR "Health Surveys"[Mesh] OR "Case Reports" [Publication Type]) AND 
(english[Language]) AND 
("2023/01/01"[Date - Publication] : "2023/02/01"[Date - Publication])
'
