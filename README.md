A STATISTICAL FRAMEWORK FOR DECOUPLING SPATIAL AND ENVIRONMENTAL DRIVERS OF 
PHENOLOGY FROM PATCHY HISTORICAL DATASETS
================================================================================

Authors: 	Swantje Grabener, Claudia Drees, Michael Staab, Wiebke Schuett, 
		Pascale Zumstein, Thorsten Assmann
Target Journal: Methods in Ecology and Evolution (Workflow Article)
Zenodo DOI:     10.5281/zenodo.XXXXXXX
License:        CC-BY 4.0 (for Data, see Acknowledgements below)


1. OVERVIEW OF THE METHODOLOGICAL FRAMEWORK
--------------------------------------------------------------------------------
Long-term ecological monitoring schemes provide vital baselines for quantifying 
climate-induced phenological shifts, yet they frequently suffer from 
discontinuous sampling effort, systematic seasonal surveillance gaps (e.g., 
winter trapping gaps), and spatiotemporal confounding between static 
geography (photoperiod and local adaptation) and dynamic meteorological drivers.

This repository provides a step-by-step Generalized Additive Mixed Modeling 
(mgcv::bam / mgcv::gam) pipeline in R that resolves these challenges across 
two core analytical stages:

  i) Effort-Standardized Baseline Template & Weighted Continuous Imputation:
     Constructs a whole-year cyclic baseline template (bs = "cc") from 
     continuous reference sequences (either internal to the monitoring network 
     or sourced from external historical baselines), derives an empirical 
     detection threshold (tau_det) to truncate inactive off-season tails to 
     structural zeros, and projects observed active-season catch rates onto the 
     baseline via ordinary least-squares through the origin (kappa). Down-
     weighted (w = 0.5) continuous pseudo-observations are injected across 
     unmonitored seasonal gaps, leveraging the non-negative real support of the 
     compound Poisson-Gamma (Tweedie) distribution to regularize cyclic splines 
     without artificial integer rounding.

 ii) Functional ANOVA Tensor Decomposition with Dynamic Environmental Drivers:
     Matches each observation day (d = DOY) with causally aligned trailing 
     meteorological drivers standardized against the network-wide daily 
     climatology, W*(l, y, d). Using orthogonal tensor-product interactions 
     (ti), the pipeline decouples geographic/photoperiodic shape modifications 
     ti(DOY, Space) from intra- and interannual meteorological plasticity 
     ti(DOY, W*) while preserving strictly isothermal comparisons across space. 
     Predictions are normalized onto the unit simplex via unconditional 
     Bayesian posterior simulation (R = 1,000 draws) to isolate phenological 
     timing from population abundance fluctuations.


2. REPOSITORY CONTENTS
--------------------------------------------------------------------------------
The repository is streamlined into a single data archive and a single 
self-contained workflow script:

  - Phenology_Pipeline_Data.RData    : Compressed R data archive (all inputs)
  - Phenology_Decoupling_Pipeline.R  : Annotated end-to-end R workflow script
  - README.txt                       : Repository documentation and usage guide


3. DESCRIPTION OF THE DATA ARCHIVE (Phenology_Pipeline_Data.RData)
--------------------------------------------------------------------------------
Loading "Phenology_Pipeline_Data.RData" in R populates four objects in the 
global environment:

  1) ecn_species_data
     A named list of three data frames ($Pterostichus_madidus, 
     $Notiophilus_biguttatus, and $Calathus_micropterus) extracted from the UK 
     Environmental Change Network (ECN) pitfall trapping dataset (1992-2015; 
     Rennie et al., 2015), including all true zero catches:
       - Site      : ECN monitoring site name (12 sites across the UK).
       - SubSite   : Replicated trapping transect identifier (47 subsites).
       - OpenDate  : Trap deployment start date (YYYY-MM-DD).
       - CloseDate : Trap collection end date (YYYY-MM-DD).
       - NTraps    : Number of active pitfall traps during the sampling interval.
       - N         : Physical individual count recorded for the focal species.

  2) wijster_nb
     Historical continuous year-round pitfall trapping records for Notiophilus 
     biguttatus from the Biological Station Wijster (1959-1966; Gerrits et al., 
     2025), used to construct the external biological baseline template:
       - StartDate : Trap deployment start date (YYYY-MM-DD).
       - EndDate   : Trap collection end date (YYYY-MM-DD).
       - N         : Physical individual count recorded.

  3) therm_daily
     Pre-calculated, causally aligned dynamic thermal drivers for every 
     (Site, Year, DOY) combination, standardized by the network-wide daily 
     climatological mean and standard deviation (derived from UK Met Office 
     daily surface air temperature records; raw temperature measurements are 
     omitted in accordance with data licensing terms):
       - Site        : ECN monitoring site name.
       - Year        : Calendar year.
       - DOY         : Day of the Year (1 to 365).
       - T_roll30    : Standardized trailing 30-day mean daily air temperature 
                       over days [d-29, d] (proximal acute thermal window).
       - GDD5_roll60 : Standardized trailing 60-day cumulative growing degree-
                       days (> 5 °C) over days [d-59, d] (developmental window).
       - GDD5_YTD    : Standardized year-to-date cumulative growing degree-days 
                       (> 5 °C) from Jan 1 to day d (ontogenetic thermal sum).

  4) site_coords
     Geographic latitude (Latitude, decimal degrees) for each of the 12 ECN 
     monitoring sites.


4. SYSTEM REQUIREMENTS & QUICK START
--------------------------------------------------------------------------------
Software & Packages:
The workflow was developed and tested in R (version 4.6.1) and requires the 
following CRAN packages:

    install.packages(c("dplyr", "lubridate", "ggplot2", "mgcv", "patchwork"))

Running the Pipeline:
Set your R working directory to the repository folder and run:

    source("Phenology_Decoupling_Pipeline.R")

Expected Runtime:
Approximately 15 to 30 seconds in total for all three species (9 candidate 
GAMMs + 3,000 posterior simulations + figure rendering) on standard desktop 
hardware, enabled by fast Restricted Maximum Likelihood (method = "fREML") 
and discrete covariate binning (discrete = TRUE) in mgcv::bam().

Outputs Produced in Working Directory:
  1. Table_1_Model_Selection_Summary.csv          (Reproduces Table 1)
  2. Champion_GAMM_PM.rds, _NB.rds, _CM.rds       (Fitted GAMM objects)
  3. Figure_2_Pterostichus_madidus_Plasticity.png (Reproduces Figure 2)
  4. Figure_4_Notiophilus_biguttatus_Plasticity.png (Reproduces Figure 4)
  5. Figure_5_Calathus_micropterus_Plasticity.png (Reproduces Figure 5)


5. ADAPTING THE WORKFLOW TO CUSTOM DATASETS
--------------------------------------------------------------------------------
To apply run_phenology_workflow() to other taxonomic groups or monitoring 
schemes (e.g., herbarium records, bird ringing campaigns, pollinator transects, 
or plankton time series):

  1. Format Observation Records: Provide a data frame containing location IDs 
     (Site, SubSite), sampling window dates (OpenDate, CloseDate), sampling 
     effort (Effort, e.g., trap-days, transect length, or observer hours), and 
     observed counts (N, including zero-count surveys).
  2. Supply a Baseline Reference: If your primary dataset contains at least one 
     continuous year-round monitoring sequence (N >= 10), the pipeline will 
     automatically detect and use it as the Step 1 baseline template. Otherwise, 
     pass an external reference dataset via the external_df argument (or inject 
     binary dormancy zeroes during known biological inactivity windows).
  3. Compute Dynamic Environmental Drivers: Calculate trailing environmental 
     metrics relevant to your study organism up to each observation day d, and 
     standardize each metric by its network-wide mean and standard deviation on 
     calendar day d.
  4. Select Model Modularity (model_type): Use model_type = "dual" for datasets 
     spanning broad geographic gradients, or model_type = "single" for single-
     site or geographically restricted studies (Manuscript Section 2.5).


6. ACKNOWLEDGEMENTS, DATA SOURCES & LICENSING
--------------------------------------------------------------------------------
  - UK Environmental Change Network (ECN): Ground beetle pitfall trapping data 
    were provided by the UK ECN and are distributed under the UK Open Government 
    Licence / NERC Environmental Information Data Centre terms:
    Rennie, S., Adamson, J., Anderson, R., Andrews, C., Bater, J., Bayfield, N. 
    et al. (2015). UK Environmental Change Network (ECN) carabid beetle data: 
    1992-2015. NERC Environmental Information Data Centre. 
    https://doi.org/10.5285/8385f864-dd41-410f-b248-028f923cb281

  - Biological Station Wijster: Historical continuous carabid trapping data 
    are licensed under CC-BY 4.0 and sourced from:
    Gerrits, G.M., Vermeulen, R., & van Klink, R. (2025). Sample event data on 
    ground beetles (Coleoptera, Carabidae) collected by Biological Station 
    Wijster (BSW) in the years 1959, 1961, 1963, 1965 and 1966. Biodiversity 
    Data Journal, 13, e171473. https://doi.org/10.3897/BDJ.13.e171473

  - UK Met Office: Daily temperature data used to derive the standardized 
    thermal indices were provided by the Met Office, UK (© Crown Copyright, 
    National Meteorological Library and Archive).

  - Funding: This work was supported by the Open Access Publication Fund of 
    Leuphana University Lüneburg.
================================================================================
