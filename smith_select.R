### Preparation

#load/install a packages
source("install_R_packages.R")

#create directories
if(!dir.exists("Ergebnisse")){dir.create("Ergebnisse")}
if(!dir.exists("errors")){dir.create("errors")}
if(!dir.exists("Bundles")){dir.create("Bundles")}


#source config
if(file.exists("config.R")&&!dir.exists("config.R")){
  source("config.R")
}else{
  source("config.R.default")  
}

#remove trailing slashes from endpoint
base <- if(grepl("/$", base)){strtrim(base, width = nchar(base)-1)}else{base}


brackets = c("[", "]")
sep = " || "

###Get all Observations between 2019-01-01 and 2021-12-31 with loinc 33763-4,71425-3,33762-6,83107-3 or 83108-1
#also get associated patient resources
#Observations have to implement MII profile
obs_request <- fhir_url(url = base, 
                        resource = "Observation", 
                        parameters = c("code" = "http://loinc.org|33763-4,http://loinc.org|71425-3,http://loinc.org|33762-6,http://loinc.org|83107-3,http://loinc.org|83108-1",
                                       "date" = "ge2019-01-01",
                                       "date" = "le2021-12-31",
                                       "_include" = "Observation:patient",
                                       "_profile" = "https://www.medizininformatik-initiative.de/fhir/core/modul-labor/StructureDefinition/ObservationLab"))

#download bundles
obs_bundles <- fhir_search(request = obs_request,
                           log_errors = "errors/observation_error.xml")

#save for checking purposes
fhir_save(bundles = obs_bundles, directory = "Bundles/Observations")

#flatten
obs_description <- fhir_table_description("Observation", 
                                          cols = c(NTproBNP.date = "effectiveDateTime",
                                                   subject = "subject/reference",
                                                   NTproBNP.value = "valueQuantity/value",
                                                   NTproBNP.unit = "valueQuantity/code",
                                                   NTproBNP.unitSystem = "valueQuantity/system"))

pat_description <- fhir_table_description("Patient",
                                          cols = c(id = "id",
                                                   gender = "gender", 
                                                   birthdate = "birthDate"))

obs_tables <- fhir_crack(obs_bundles, 
                         design = fhir_design(obs = obs_description, pat = pat_description),
                         data.table = TRUE)

if(nrow(obs_tables$obs)==0){
  write("Konnte keine Observations für NTproBNP auf dem Server finden. Abfrage abgebrochen.", file ="errors/error_message.txt")
  stop("No NTproBNP Observations found - aborting.")
}

if(nrow(obs_tables$pat)==0){
  write("Konnte keine Patientenressourcen für NTproBNP-Observations auf dem Server finden. Abfrage abgebrochen.", file ="errors/error_message.txt")
  stop("No Patients for NTproBNP Observations found - aborting.")
}

### get associated resources: consent, encounters, conditions

#split patient id list into smaller chunks that can be used in a GET url 
#(split because we don't want to exceed allowed URL length)
patients <- obs_tables$pat$id #all patient ids
nchar_for_ids <- 1800 - nchar(paste0(base, "Encounter?_profile=https://www.medizininformatik-initiative.de/fhir/core/modul-fall/StructureDefinition/KontaktGesundheitseinrichtung")) #assume maximal length of 1800

n <- length(patients)
list <- split(patients, ceiling(seq_along(patients)/n)) #split patients ids in chunks of size n
nchar <- sapply(list, function(x){sum(nchar(x))+(length(x)-1)}) #compute number of characters for each chunk, including commas for seperation

#reduce the chunk size until number of characters is small enough
while(any(nchar > nchar_for_ids)){
  n <- n/2
  list <- split(patients, ceiling(seq_along(patients)/n))
  nchar <- sapply(list, function(x){sum(nchar(x))+(length(x)-1)})
}


#get consent
if(filterConsent){

  consent_list <- lapply(list, function(x){

    ids <- paste(x, collapse = ",")

    consent_request <- fhir_url(url = base,
                            resource = "Consent",
                            parameters = c(patient = ids))

    consent_bundles <- fhir_search(consent_request,
                               log_errors = "errors/consent_error.xml")

  })
  #bring consent results together, save and flatten
  consent_bundles <- fhircrackr:::fhir_bundle_list(unlist(consent_list, recursive = F))
  fhir_save(bundles = consent_bundles, directory = "Bundles/Consents")
  
  consent_description <- fhir_table_description("Consent",
                                                cols = c(patient = "patient/reference",
                                                         provision.code = "provision/provision/code/coding/code",
                                                         provision.system = "provision/provision/code/coding/system"))
  consent_table <- fhir_crack(consent_bundles, 
                              design = consent_description,
                              data.table = TRUE,
                              brackets = brackets,
                              sep = sep)
  
  #unpack multiple provision info
  consent_table <- fhir_melt(consent_table,columns = c("provision.code", "provision.system"), brackets = brackets, sep = sep, all_columns = T)
  consent_table <- fhir_melt(consent_table,columns = c("provision.code", "provision.system"), brackets = brackets, sep = sep, all_columns = T)
  consent_table <- fhir_rm_indices(consent_table, brackets = brackets)
  
  #find Patients that have code for MDAT_wissenschaftlich_nutzen_EU_DSGVO_NIVEAU
  allowed_pats <- consent_table[provision.code=="2.16.840.1.113883.3.1937.777.24.5.3.8" & provision.system=="urn:oid:2.16.840.1.113883.3.1937.777.24.5.3"]
  allowed_pats <- sub("Patient/", "", allowed_pats$patient)
  
}

###merge observation and patient data
#prepare key variables for merge
obs_tables$obs[, subject:=sub("Patient/", "", subject)] 

#sort out col types
obs_tables$obs[, NTproBNP.date := as.Date(NTproBNP.date)]

#merge
obsdata <- merge.data.table(x = obs_tables$obs, 
                            y = obs_tables$pat, 
                            by.x = "subject",
                            by.y = "id",
                            all.x = TRUE)


#if necessary filter for consent and create new list of patient id chunks
if(filterConsent){
  obsdata <- obsdata[subject %in% allowed_pats]
  
  #split patient id list into smaller chunks that can be used in a GET url 
  #(split because we don't want to exceed allowed URL length)
  patients <- obsdata$subject #filtered patient ids
  nchar_for_ids <- 1800 - nchar(paste0(base, "Encounter?_profile=https://www.medizininformatik-initiative.de/fhir/core/modul-fall/StructureDefinition/KontaktGesundheitseinrichtung")) #assume maximal length of 1800
  
  n <- length(patients)
  list <- split(patients, ceiling(seq_along(patients)/n)) #split patients ids in chunks of size n
  nchar <- sapply(list, function(x){sum(nchar(x))+(length(x)-1)}) #compute number of characters for each chunk, including commas for seperation
  
  #reduce the chunk size until number of characters is small enough
  while(any(nchar > nchar_for_ids)){
    n <- n/2
    list <- split(patients, ceiling(seq_along(patients)/n))
    nchar <- sapply(list, function(x){sum(nchar(x))+(length(x)-1)})
  }
  
}

#get encounters
encounter_list <- lapply(list, function(x){
  
  ids <- paste(x, collapse = ",")
  
  enc_request <- fhir_url(url = base,
                          resource = "Encounter",
                          parameters = c(subject = ids,
                                         "_profile" = "https://www.medizininformatik-initiative.de/fhir/core/modul-fall/StructureDefinition/KontaktGesundheitseinrichtung"))

  enc_bundles <- fhir_search(enc_request,
                             log_errors = "errors/encounter_error.xml")

})

#bring encounter results together, save and flatten
encounter_bundles <- fhircrackr:::fhir_bundle_list(unlist(encounter_list, recursive = F))
fhir_save(bundles = encounter_bundles, directory = "Bundles/Encounters")

enc_description <- fhir_table_description("Encounter",
                                          cols = c(subject = "subject/reference",
                                                   encounter.start = "period/start", 
                                                   encounter.end = "period/end",
                                                   serviceType = "serviceType"))
enc_table <- fhir_crack(encounter_bundles, 
                        design = enc_description,
                        data.table = TRUE)
if(nrow(enc_table)==0){
  write("Konnte keine Encounter-Ressourcen zu den gefundenen Patients finden. Abfrage abgebrochen.", file ="errors/error_message.txt")
  stop("No Encounters for Patients found - aborting.")
} 


#get conditions
condition_list <- lapply(list, function(x){
  
  ids <- paste(x, collapse = ",")

  con_request <- fhir_url(url = base,
                          resource = "Condition",
                          parameters = c(subject = ids,
                                         "_profile" = "https://www.medizininformatik-initiative.de/fhir/core/modul-diagnose/StructureDefinition/Diagnose"))
  
  con_bundles <- fhir_search(con_request,
                             log_errors = "errors/condition_error.xml")
  
})

#bring together conditions, save and flatten
condition_bundles <- fhircrackr:::fhir_bundle_list(unlist(condition_list, recursive = F))
fhir_save(bundles = condition_bundles, directory = "Bundles/Conditions")

con_description <- fhir_table_description("Condition",
                                          cols = c(clinicalStatus.code = "clinicalStatus/coding/code",
                                                   clinicalStatus.system = "clinicalStatus/coding/system",
                                                   verificationStatus.code = "verificationStatus/coding/code",
                                                   verificationStatus.system = "verificationStatus/coding/system",
                                                   code = "code/coding/code",
                                                   code.system = "code/coding/system",
                                                   subject = "subject/reference",
                                                   onsetPeriod.start = "onsetPeriod/start",
                                                   onsetPeriod.end = "onsetPeriod/end",
                                                   recordedDate = "recordedDate"))

con_table <- fhir_crack(condition_bundles, 
                        design = con_description,
                        data.table = TRUE,
                        brackets = brackets,
                        sep = sep)



###merge encounter data
#prepare key variables for merge
enc_table[, subject:=sub("Patient/", "", subject)]

#sort out col types
enc_table[, encounter.start := as.Date(encounter.start)]
enc_table[, encounter.end := as.Date(encounter.end)]

#merge
cohort <-  merge.data.table(x = enc_table, 
                                   y = obsdata, 
                                   by.x = "subject",
                                   by.y = "subject",
                                   all.x = TRUE)

#filter encounters that dont belong to the NTproBNP observation
cohort <- cohort[NTproBNP.date >= encounter.start & NTproBNP.date <= encounter.end]

if(nrow(con_table) > 0){
  #filter conditions
  #expand codes
  conditions <- fhir_melt(con_table, columns = c("code", "code.system"), brackets = brackets, sep = sep, all_columns = TRUE)
  conditions <- fhir_melt(conditions, columns = c("code", "code.system"), brackets = brackets, sep = sep, all_columns = TRUE)
  conditions <- fhir_rm_indices(conditions, brackets = brackets)
  conditions[,resource_identifier:=NULL]
  
  #filter for ICD codesystem
  conditions <- conditions[grepl("icd-10", code.system)]
  
  #codes for selection
  condition_codes <- c("I48.0", "I48.1", "I48.2", "I48.9", "I50.00", "I50.01", "I50.02",
                       "I50.03", "I50.04", "I50.05", "I50.1", "I50.11", "I50.12", "I50.13",
                       "I50.14", "I50.19", "I50.9", "I13.21", "I13.20", "I13.3", "I13.01", "I13.00")
  
  #filter for specific codes
  conditions <- conditions[code %in% condition_codes]
  
  #prepare key variable
  conditions[, subject:=sub("Patient/", "", subject)]
}else{
  conditions <- con_table
}


###Export
if(!dir.exists("Ergebnisse")){dir.create("Ergebnisse")}
write.csv2(cohort, paste0("Ergebnisse/Kohorte.csv"))
write.csv2(conditions, paste0("Ergebnisse/Diagnosen.csv"))





