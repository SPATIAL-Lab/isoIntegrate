# isoIntegrate
Evaluation of residential mobility and forward modeling of oxygen isotope integration in human enamel. Here, we evaluate how childhood residential mobility affects geographic provenance interpretations based on stable oxygen isotopes in human enamel and present first-order, residence-informed forward models for predicting the integrated isotope signal recovered from conventional forensic samples.

# code/
Contains all code to reproduce results presented in Delgado et al. in prep.
## X_setup.R
Reads in required data and loads functions required for data analysis. Functions built without library dependencies to ensure stability throughout the analysis. Functions include:
* _rmse_- Calculation for root mean squared error
* _mae_ - Calculation for mean absolute error
* _bootstrap_error_single_ - Single dataset bootstrap for error metric calculation and confidence interval determination
* _bootstrap_error_pair_ - Pairwise bootstrap for comparison of error metrics between forward model residuals
* _make_boot_df_ - Cleaning function for paired bootstrap output to manipulate data for visualization
* _glmm_summary_ - Creates a coefficient table from GLMM outputs

## Y_analysis.R
Performs all data analysis as presented in Delgado et al. in prep. Sources X_setup. This includes: 
* Regression and predictions of isotopic data.
* Forward modeling of integrated isotope measurements from residence histories
* Diagnostic statistics (Shapiro-Wilk Normality, Levene's Test).
* Error metric calculation and bootstrap for individual datasets and pairwise comparisons between forward models.
* Generalized linear mixed effects models for absolute prediction error as a response to residence history attributes (including one sensitivity test) and code supporting their interpretation and pairwise comparisons of estimated marginal means.
Option to load in exact coefficient table used in initial interpretations; however, the GLMM outputs will be stable across iterations of the script but this allows a bypass for running the entire script if desired.

## Z_figures.R
Recreates the figures present in Delgado et al. in prep. Sources Y_analysis.

Option to load in exact pairwise bootstrap results to visualize these comparisons; however, seeds within Y_analysis should be stable for reproduction.

# data/
Contains the data required for the analysis. Includes:
* nonmovers_deident.csv - Deidentified FIND-EM data for nonmobile participants. **participant_id** defines the unique participant identified; **n.res** defines the number of residences within the span of the modeled crown formation interval; **res_type** defines the type of residence over the span of survey data; **timelength_beg** defines the age (years) in which the residence began and is constrained to the modeled crown formation interval; **timlength_end** defines the age (years) in which the residence ends and is constrained to the modeled crown formation interval; **state** defines the state of residence; **timelength** defines the total number of years of the formation interval resided at the residence. **d18O.m** defines the mean measured oxygen isotope ratio within the participant's enamel carbonate; **d18O.sd** defines the standard deviation of the measured oxygen isotope ratio within the participant's enamel carbonate.
* movers_deident.csv - Deidentified FIND-EM data for mobile participants. Includes two additional variables: **avg_distance** which defines the average distance (miles) between the residences and **d18O.twzcta.mean** which defines the expected mean tap water oxygen isotope composition for a given ZCTA as modeled in the isoscape from Verostick et al. (2025).
* dest_orig.csv - Cleaned and manipulated census data on 2023 residential mobility. Used to create the migration_map figure and not core for analysis.
* enamLP.png - Enamel carbonate d18O isoscape linear regression predictions informed by tap water isoscape. Used to create the migration_map figure and not core for analysis.
* coef_big.rds - R data format for the GLMM coefficient table generated in the initial analysis. For reference.
* pairwise_results.rds - R data format for the pairwise bootstrap error metric comparison for forward models. For reference.
* states.shp/.shx/.dbf - Shapefile for CONUS state boundaries used to create the migration_map figure and not core for analysis.

# Figures/
Contains the figures presented in Delgado et al. in prep. Reference manuscript for figure descriptions.

# Tables/
Contains the formatted Excel tables presented in Delgado et al. in prep including error metrics, GLMM coefficient table, GLMM metadata, and sensitivity test comparison.

