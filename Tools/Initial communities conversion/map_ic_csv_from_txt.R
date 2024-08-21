library("stringr")

#define some functions
fill_NA <- function(x) {
  which.na <- c(which(!is.na(x)), length(x) + 1)
  values <- na.omit(x)
  
  if (which.na[1] != 1) {
    which.na <- c(1, which.na)
    values <- c(values[1], values)
  }
  
  diffs <- diff(which.na)
  return(rep(values, times = diffs))
}

get_biomass <- function(biomass_vec){
  biomass <- biomass_vec[1, ]%>%
    gsub(" ", "", .) %>%
    gsub("\\(", "", .) %>%
    gsub("\\)", "", .)
}

get_age <- function(age_vec){
  age <- age_vec[1, ]%>%
    gsub(" ", "", .) %>%
    gsub("\\(", "", .) %>%
    gsub("\\)", "", .)
}

get <- `[` #dumb workaround for using the built-in pipes with special functions


#read and clean up data
ic <- readLines("initial-communities_LTB.txt")

ic <- gsub("\\s+", " ", ic) #replace all whitespace with individual space
ic <- gsub("^\\s+", "", ic) #remove whitespace from beginning of strings
ic <- ic[!(ic == "")]
ic <- ic[!(grepl(">>", ic))]
ic <- ic[!(grepl("LandisData", ic))]

mapcode_lines <- as.integer(sub(".*?MapCode.*?(\\d+).*", "\\1", ic)) |> fill_NA()

ic_frame <- data.frame(mapcode_lines, ic) |> get(!grepl("MapCode", ic), )
ic_frame$SpeciesCode = stringr::str_extract(ic_frame$ic, '\\w*')

#get biomass and ages
biomass_regex <- regmatches(ic_frame$ic, gregexec(pattern = "(\\(\\d+\\))", ic_frame$ic)) %>%
  lapply(get_biomass)
age_regex <- regmatches(ic_frame$ic, gregexec(pattern = "(\\d+ \\()", ic_frame$ic)) %>%
  lapply(get_age)

lengths <- lapply(biomass_regex, length)  #how many cohorts per species/site

ic_new <- data.frame(MapCode = rep(ic_frame$mapcode_lines, times = lengths),
                     SpeciesName = rep(ic_frame$SpeciesCode, times = lengths),
                     CohortAge = unlist(age_regex),
                     CohortBiomass = unlist(biomass_regex)
                     )
write.csv(ic_new, "ic.csv")


