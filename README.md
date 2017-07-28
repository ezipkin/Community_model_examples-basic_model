# Community model examples
# Part I - Basic model
Based on work from *Zipkin E.F., Royle J.A., Dawson D.K., and Bates S. 2010. Multi-species occurrence models to evaluate the effects of conservation and management actions. Biological Conservation. 143: 479-484.*

## **Description:**
The hierarchical community model is a multi-species approach to obtain community information, such as species or assemblage richness, by estimating individual species occurrence probabilities. The fundamental idea behind the community approach is that collective information on all observed species can inform probabilities of detection and occurrence for both observed and unobserved species, even those that are rare or elusive. This results in an improved composite analysis of the community and increased precision in species specific estimates of occurrence. The hierarchical model can be specified to incorporate habitat and sampling effects that influence occurrence and detection. Thus the community approach can provide the best possible estimates of species richness and other metrics of interest across a heterogeneous landscape, while accounting for variation in occurrence and detection among species.

This repo provides R and WinBUGS code to implement a "basic" version of this modeling framework. The code here estimates species-specific detection and occupancy assuming no differences among locations. It also estimates species richness using data augmentation. The model assumes no covariates on either the biological (occupancy) or sampling processes.

See: https://www.mbr-pwrc.usgs.gov/site/communitymodeling/home/ for more information.

## **Data:**
occ data.csv - Bird species detection nondetection data. Each column is an observation containing information on 1) the site where the observation occured (Point), 2) the time (Time) and date (Date) of the detection, 3) the species that was detected (Species), and 4) which survey replicate (Rep) the observation data was collected (as determined by the unique site/date combination).


## **Code:**
basic model code.R - R code to run the multi-species occupancy model with no covariates (uses the data in the occ_data.csv file). Contains code to import and reshape the data, create the BUG model file, and run the model file in WinBUGS. There is also code to process the results and make some figures.
