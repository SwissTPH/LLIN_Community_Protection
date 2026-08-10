# Load the necessary packages
library(devtools)
library(openMalariaUtilities)
library(AnophelesModel)
library(tidyverse)
library(omucompat)

#########################################################3
git_repo_path="/scicore/home/pothin/suboi0000/GitRepositories/om-ghana"
setwd(file.path(git_repo_path,"emulator_training/OMsimulations/nsp2023/COMMUNITY_PROTECTION/COHORTS/NEW_LLIN_PARAMETERIZATION/Final_Simulation"))

ITNcascades_path="/scicore/home/pothin/GROUP/OM_countrymodelling/Helpers/ITNcascades-main/RCT_validation"
vectorParams = "/scicore/home/pothin/GROUP/OM_countrymodelling/Helpers/VectorControl_template"

parameters_GVI=read.csv(file.path(ITNcascades_path, "../fitted_parameters_posteriormax_ITN_FULL.csv")) #|> as.data.frame()
source(file.path(vectorParams,"helpfunction_omu_gvi.R"))
source(file.path(git_repo_path, "emulator_training/OMsimulations/nsp2023/COMMUNITY_PROTECTION/COHORTS/NEW_LLIN_PARAMETERIZATION/get_activity_per_species.R") )



# Function that creates a list with all the elements which are specific for the
# base xml (placeholders, interventions, etc.)
create_baseList = function(country_name="GHA", 
                           sim_start="1918-01-01", 
                           versionnum=44L) {
  
  ## Basic xml skeleton
  baseList = list(
    # Mandatory
    expName = country_name,
    # Mandatory
    OMVersion = versionnum,
    # Mandatory
    demography = list(),
    monitoring = list(),
    interventions = list(),
    healthSystem = list(),
    entomology = list(),
    # These are optional for OM
    # parasiteGenetics = list(),
    # pharmacology = list(),
    # diagnostics = list(),
    model = list()
  )
  
  baseList = defineDemography(
    baseList,
    name = country_name,
    popSize = "@pop@",
    maximumAgeYrs = 94,
    lowerbound = 0,
    poppercent = GHA$poppercent,
    upperbound = GHA$upperbound
  )
  
  ## Create monitoring snippet
  baseList[["monitoring"]] = list(
    name = "Surveys",
    ## Mandatory, different from OM schema
    startDate = sim_start,
    continuous = monitoringContinuousGen(period = 1,
                                         options = list(
                                           name = c("human infectiousness", "N_v0",#"input EIR", "simulated EIR",
                                                    "immunity h", "immunity Y", 
                                                    "new infections",
                                                    "num transmitting humans", 
                                                    "ITN coverage", "GVI coverage", 
                                                    "alpha", "P_B", "P_C*P_D"),
                                           value = c( "true", 
                                                     "true", "true", "true", 
                                                     "true", "true",
                                                     "true", "true", 
                                                     "true", "false", "false") #"true", "true",
                                         )
    ),
    SurveyOptions = monitoringSurveyOptionsGen(
      options = list(
        name = c("nHost", "nPatent", "nUncomp", "nSevere", "nDirDeaths", 
                 "nInfect", "innoculationsPerAgeGroup",
                 "nTreatments1","nTreatments2",
                 "nTreatments3","expectedDirectDeaths"), #"inputEIR", "simulatedEIR",
        value = c("true", "true", "true", "true", "true", "true", "true", 
                  "true", "true", "true", "true")# , "true", "true"
      )
    ),
    # We are setting a survey on the 5th day of each month
    surveys = monitoringSurveyTimesGen(detectionLimit = 100, #RDT detection limit
                                       startDate = "2015-01-01", #monitoring start date - 2015
                                       endDate = "2035-01-01",
                                       interval = list(days = c(5), 
                                                       months = c(1:12), 
                                                       years = c(2015:2035)),
                                       simStart = sim_start),
    ## surveyAgeGroupsGen will write thirdDimension table to cache, important for postprocessing
    ageGroup = surveyAgeGroupsGen(lowerbound = 0, upperbounds = c(1, 2, 5, 10, 100)),
    # add cohorts in the monitoring snippet
    cohorts = list(subPop = list(id="LLINusers", number="1"))
  )
  
  # which mosquitoes will be used ?
  mosqs  = c("funestus_indoor1","funestus_indoor2","funestus_outdoor",
             "gambiae_indoor1","gambiae_indoor2","gambiae_outdoor",
             "arabiensis_indoor1","arabiensis_indoor2","arabiensis_outdoor" )
  
  # contrib: the relative contribution of each 'mosquito' to EIR
  contrib = c("@fin1@","@fin2@","@fout@",
              "@gin1@","@gin2@","@gout@",
              "@ain1@","@ain2@","@aout@")
  
  cbind( mosqs, contrib)
  ##############################################################
  #############
  # RELATIVE CONTRIBUTION OF DIFFERENT SPECIES for GH
  
  #NMESP - 95% gambiae and funestus
  #Colemane et al: https://link.springer.com/article/10.1186/s13071-023-05793-2
  #props - 83% gambiae
  #implies 12% funestus
  #assume 5% arabiensis or remove arabiensis
  
  #sporozoite rates
  #Osae et al. - https://pmc.ncbi.nlm.nih.gov/articles/PMC11107868/?utm_source=chatgpt.com#ref1
  #dominant species are funestus and gambiae
  #funestus=0.049, gambiae=0.038
  
  # function to reweight species proportions with spoporzoite rate.
  # prop_raw: raw proportions of each vector species
  # sporozoite_rates: sporozoite rates for each species. Default values are from Kweyamba et al. 2025 https://www.nature.com/articles/s41598-025-86409-w
  
  # sporozoite_rate_funestus=0.049
  # sporozoite_rate_gambiae=0.038
  # sporozoite_rate_arabiensis=0.0
  
  #use default since this is a more generic experiment
  sporozoite_rate_funestus=0.235
  sporozoite_rate_gambiae=0.114
  sporozoite_rate_arabiensis=0.049
  
  reweight_species_contributions=function(prop_raw=c("prop_funestus"=NA, "prop_gambiae"=NA, "prop_arabiensis"=NA),
                                          sporozoite_rates=c("sporozoite_rate_funestus"=0.235, "sporozoite_rate_gambiae"=0.114, "sporozoite_rate_arabiensis"=0.049)){
    
    
    if(sum(prop_raw)!=1){
      stop("the raw proportions don't sum up to 1")
    }
    contri_funestus=prop_raw["prop_funestus"]*sporozoite_rates["sporozoite_rate_funestus"]
    contri_gambiae=prop_raw["prop_gambiae"]*sporozoite_rates["sporozoite_rate_gambiae"]
    contri_arabiensis=prop_raw["prop_arabiensis"]*sporozoite_rates["sporozoite_rate_arabiensis"]
    
    tot=contri_funestus+contri_gambiae+contri_arabiensis
    return(c("prop_funestus"=as.numeric(contri_funestus/tot), "prop_gambiae"=as.numeric(contri_gambiae/tot), "prop_arabiensis"=as.numeric(contri_arabiensis/tot)))
    
  }
  
  
  # enter here the initial proportions for each species
  prop_update=reweight_species_contributions(prop_raw=c("prop_funestus"=0.3, "prop_gambiae"=0.6, "prop_arabiensis"=0.1))
  
  ###############
  # INDOOR/ OUTDOOR
  
  # calculate vector susceptibility to intervention
  
  all_activity=split_complex_mean(activity_pat_new)
  activity_gambiae=all_activity$gambiae_ss
  activity_funestus=all_activity$funestus
  activity_arabiensis=all_activity$funestus
  
  get_indoor1_indoor2=function(activity, mosquito_species){
    ent_params= def_vector_params(mosquito_species = mosquito_species)
    exposure_vector=get_in_out_exp(vec_p =ent_params,  activity_cycles = activity )
    
    indoor1=min(exposure_vector$Exposure_Indoor_whileinbed, exposure_vector$indoor_resting)
    indoor2=max(exposure_vector$Exposure_Indoor_whileinbed, exposure_vector$indoor_resting)-min(exposure_vector$Exposure_Indoor_whileinbed, exposure_vector$indoor_resting)
    outdoor=1-indoor2-indoor1
    
    which_min=ifelse(exposure_vector$Exposure_Indoor_whileinbed<= exposure_vector$indoor_resting,
                     "in bed", "indoor biting")
    
    return(list("exp_vector"=c("indoor1"=indoor1, "indoor2"=indoor2, "outdoor"=outdoor),
                which_is_indoor1=which_min))
  }
  
  
  exposure_gambiae=get_indoor1_indoor2(activity=activity_gambiae, mosquito_species = "Anopheles gambiae")
  exposure_funestus=get_indoor1_indoor2(activity=activity_funestus, mosquito_species = "Anopheles funestus")
  exposure_arabiensis=get_indoor1_indoor2(activity=activity_arabiensis, mosquito_species = "Anopheles arabiensis")
  
  ##############################################################
  
  ## Entomology section
  baseList = make_ento_compat(baseList = baseList, mosqs, contrib, EIR = "@EIR@",
                              seasonality = paste0("@m", 1:12, "@"), ## NEED TO MATCH with
                              propInfected = .078, propInfectious = .021)
  
  ## Specify seed and finish XML file
  baseList = write_end_compat(baseList = baseList,
                              seed = "@seed@", modelname = "base")  
  
  ### ALL HUMAN INTERVENTIONS NEED TO BE REVIEWED
  # Begin interventions for humans
  ########################################################################################################
  ## Definition section
  ########################################################################################################
  
  # No histITNs to be deployed
  
  # 1 step: Define a dummy historical ITN in the baselist
  baseList <- define_ITN(
    baseList = baseList, component = "dummyITN", mosquitos = mosqs, historical = TRUE,
    resist = TRUE, halflife = 0, strong=FALSE
  )
  
  # 2 step: define your historical ITNs with the GVI function
  distrib_year=2019
  distrib_month=1
  # 2a) step: define the parameters for your Weibull function (half life and the shape parameter) 
  #to be used in both the historical and future GVI snippets
  
  halflife_IG2=2.4
  kappa_IG2=2.4
  
  
  # 2b) step: define the extract_GVI_params which helps you to extract parameters for the specific EHT
  ## West Africa CI EHT
  
  ## West Africa CI EHT - only IG2 deployed
  IG2_AssengaCI = extract_GVI_params(EHT="Assenga, CI",
                                     netType="Interceptor G2",
                                     halflife_functionalSurvival=halflife_IG2,
                                     kappa_functionalSurvival=kappa_IG2,
                                     myname="IG2_AssengaCI",
                                     mosqs=mosqs,
                                     parameters_GVI=parameters_GVI,
                                     insecticide_decay = F,
                                     decay = "step",
                                     active_categories = c("indoor1", "indoor2"))#
  
  
  baseList= defineGVI_simple(baseList = baseList,
                             vectorInterventionParameters=IG2_AssengaCI,
                             append = FALSE,
                             verbatim = TRUE,
                             hist = FALSE)
  
  print(paste("Finished setting up nets") )
  
  ## Define empty ITN/IRS   #### This is to have something when coverage is zero??
  baseList <- define_nothing_compat(
    baseList = baseList, component = "nothing", mosqs = mosqs)
  
  
  ########################################################################################################
  ## Deployment section
  ########################################################################################################
  
  
  ############################
  ## Future  (TO BE REVIEWED)
  ############################
  
  # from Clara, no need for the deployment below, use what is in the wiki (sent on Teams chat) 
  
  # ITN_years <- c(2023, 2026, 2029)
  list_ITN_snippets_id=c("IG2_AssengaCI")
  
  for(GVI_snippet_id in list_ITN_snippets_id){
    baseList <- deploy_it_compat(
      baseList = baseList, component = paste0(GVI_snippet_id),
      coverage = paste0("@futNetcovstart2023@"),
      ## allowing for different hist coverage levels
      byyear = FALSE,
      ##  annual deployments
      y1 = 2023, y2 = 2029 , 
      every = 3, 
      interval = "year",
      m1 = 1, m2 = 1, d1 = 5, d2 = 5,
      SIMSTART = sim_start
    )
  }
  
  ############################
  ## End of Future deployment
  ###########################
  
  #### HUMAN INTERVENTIONS STOP HERE
  
  
  # Importation
  baseList = define_importedInfections_compat(baseList = baseList, 10, time = 0)
  
  # Health system changes (we set the coverage of access to treatment, the effective coverage, each year) for past
  baseList = define_changeHS_compat(baseList = baseList, access = "histAccess",
                                    y1 = 2000, y2 = 2022,
                                    use_at_symbol = TRUE,
                                    ## Default values
                                    pSelfTreatUncomplicated = 0.01821375,
                                    pSeekOfficialCareSevere = .48,
                                    SIMSTART = sim_start)
  # Health system changes (we set the coverage of access to treatment, the effective coverage, each year) for future
  baseList = define_changeHS_compat(baseList = baseList, access = "futCMcov",
                                    y1=2023,y2=2031,
                                    use_at_symbol = TRUE,
                                    ## Default values
                                    pSelfTreatUncomplicated = 0.01821375,
                                    pSeekOfficialCareSevere = .48,
                                    SIMSTART = sim_start)
  
  
  # Write a "dummy" health system (the simulation will use the previously defined health system changes)
  baseList = write_healthsys_compat(baseList = baseList, access = 0)
  
  return(baseList)
}

