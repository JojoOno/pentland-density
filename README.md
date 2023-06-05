# pentland-density
Predicting the density of harbour seals across the Pentland Firth, based purely on North Coast tags, for use in collision risk models
Data and Justification
Density and distribution of seals in the UK is relatively well understood, especially in areas where large volumes of telemetry data exists to inform habitat preference models (Carter et al. 2022). However, current iterations of seal at-sea distribution predictions can be lacking in spatial resolution when used to determine the environmental impact of smaller scale human activities such early stage renewable energy installations and demonstration projects, port and harbour developments and aquaculture sites. Further, these predictions do not presently take temporal fluctuations into account, which in some more dynamic regions can result in under or overpredicting impacts which themselves vary over similar temporal scales. One such example of this is seal density in the Pentland Firth, Scotland, where seal distribution and behaviour is known to be heavily affected by tidal dynamics (Onoufriou et al. 2020; Onoufriou et al. 2021). This is of particular interest given the installation of tidal energy converters in the region which pose a potential collision risk and also vary in their operational state over a tidal cycle. Therefore, the work described below seeks to estimate harbour seal density within this region, taking variation across tidal cycles into account. 

Data from GPS telemetry devices deployed on seals hauling out on the north coast were leveraged for this purpose. Data collection protocols, seal capture sites and data details are provided in Onoufriou et al. (2021). High resolution location fixes were derived from two different tag types. Fastloc® GPS/GSM tags (SMRU Instrumentation) were deployed on 14 harbour seals in 2011 and 2012 and had a minimum inter-location interval of 5 minutes, and Fastloc® GPS/UHF tags (Pathtrack Ltd.) were deployed on 40 harbour seals over 4 deployments in 2016, 2017 and 2018 and had a minimum inter location interval of 3 minutes. Seal tracks were linearly interpolated to regularised 30-minute intervals to ensure model data resolution was consistent between all individuals. A 30-minute interval was chosen as this represented a compromise between data-resolution differences between tag types. Interpolated locations which fell within data gaps of >2 hours were determined to be unreliable and removed from the analyses. Further, data between May and September (1.3% of all locations) were removed to ensure that behavioural responses to the turbine array were not conflated with changes due to breeding and moulting at this time of the year. 

Statistical Analysis 
Telemetry data are, by their nature, presence only data. In order to analyse seal distribution with respect to tidal state, a use-availability design was used (Matthiopoulos et al. 2003). Harbour seals are central place foragers, ultimately constrained by distance to suitable haulouts (Bailey et al. 2014). A single accessibility polygon was created around the central haulout site along the north coast, the radius of which was based on the largest distance value of all trips across all individuals. For each observed presence point, two temporally matched, randomly placed pseudo-absences were generated within the polygon. Pseudo-absences can be considered a representative sample of points within the seals’ accessible range and provide a means of modelling the difference between the space observed to be used by the seals and the space that is approximately available to them. The number of pseudo-absences was restrained by the additional running time, but previous analyses on the same data have revealed that coefficients were stable beyond 2 pseudo-absences per presence point (Onoufriou et al. 2021). 

In accordance with the use-availability approach, each model was fitted using a binary response of 0=pseudo-absence and 1=presence.  In the model, distribution was represented using spatial smooths; a spatial smooth in this case represented a 2-dimensional smooth of latitude and longitude. Given the influence of tidal state on harbour seal behaviour at sea, a separate spatial smooth was offered for 4 tidal phases: low water (3 hour period centred around peak low water), high water (3 hour period centred around peak high water), flood tide (between the end of the 3-hour ‘low water’ and the 3-hour ‘high water’ periods) and ebb tide (between the end of the 3-hour ‘high water’ and the 3-hour ‘low water’ periods). Tidal state information was extracted from the tidal prediction software POLPRED version 2.003 (National Oceanography Centre, Liverpool, UK). Given the relatively small study area and short trip distances of the seals, tidal states (low and high waters) were predicted for a single point at the central location of the turbine array as tidal state would not be expected to vary considerably over the study site. To generate the required discrete tidal covariate, the difference between the time of each seal location and the nearest high-water time was calculated.  

 
The models were fitted within the R package MRSea (Scott-Hayward et al. 2013). This allowed the inclusion of Complex Region Spatial Smoothers (CReSS) which consider geodesic distances around the coast.  This avoids problems with typical spatial smoothing algorithms that smooth with Euclidean distance, often across impassable barriers, and produce poor estimates of animal densities in topographically complex areas. CReSS smooths were employed in combination with Spatially Adaptive Local Smoothing Algorithms (SALSA) to select for the most appropriate number and location of knots. This algorithm has recently been adapted for using survey data to investigate the effects of anthropogenic structures on marine species’ distributions and can be adapted  for studies using telemetry data to assess changes in distributions using use-availability designs (Russell et al. 2016). For each of 27 model runs (fitting between 4 and 30 knots) the position of the knots chosen by SALSA was that which minimised AIC. The final number of knots determined how variable the 2-dimensional surface was, with a larger number of knots indicating greater variability across space. In this case, 25 knots were retained indicating a large variation across space, with several knot locations apparent close to haulout sites (see figure: /results/knot_position.png).


Traditional generalised additive model (GAM) inference assumes independence between model residuals.  Fitting models to telemetry data often violates this assumption, given the likelihood of sequential data points being heavily dependent on temporally adjacent observations. Final models were therefore re-fitted using Generalised Estimating Equations (GEE) to ensure the results were robust to the presence of any residual autocorrelation. GEEs allow for entire time-series of data to be modelled in a regression analysis while accounting for residual auto-correlation without explicitly modelling it. This approach requires data to be split into discrete panels, between which independence is assumed but within which the autocorrelation is accounted for through robust, sandwich-based estimates of variance. Results from Onoufriou et al. (2021) suggested that individual seal was the most appropriate panel size for this analysis as it reduced autocorrelation to below statistically significant correlation at a lag of 1 i.e. correlation between neighbouring residuals was not deemed significant. As pseudo-absences are randomly generated, and thus will not be associated with residual autocorrelation, each pseudo-absence was included in separate panels before final models were run as per Russell et al. (2016).

Model selection was carried out using an assessment of marginal p-value for the interaction term. This was calculated using a Wald’s test from the ‘getPvalues’ function in the R package ‘MRSea’. Using a p-value significance threshold of 0.05, the tidal-interaction was determined to be significant. The exponent of the linear predictor from the final model was used to generate predictions of distribution. A prediction grid comprised of 500 m x 500 m grid-cells was constructed and predictions were generated for each tidal state. 

Results

The final model retained the interaction between the smooth of lat/lon location and tidal state, which is demonstrated in the plots below (see figure: /results/preds_all-states.png). Overall, there is a high degree of coastally focussed movement, with the highest densities predicted in close proximity to haul-out sites. As would be expected, these coastal distributions are not as pronounced during high tide, where available haul-out substrate is decreased, resulting in a larger degree of offshore activity. Offshore activity appears focussed to regularly used corridors between hot-spots. Aside from coastal hot-spots around haul-out sites, the predominant hot-spots exist in the far western and south-western extremities of the study area during high and ebbing tide, with a return to sites in the inner and north-eastern regions during low and flooding tides. This may be reflective of more passive movement with the prevailing tidal flow as the tide ebbs to the west (where densities are highest during high and ebbing tide). 

References

Bailey, H., P.S. Hammond, and P.M. Thompson, Modelling harbour seal habitat by combining data from multiple tracking systems. Journal of Experimental Marine Biology and Ecology, 2014. 450: p. 30-39.

Carter, M.I., Boehme, L., Cronin, M.A., Duck, C.D., Grecian, W.J., Hastie, G.D., Jessopp, M., Matthiopoulos, J., McConnell, B.J., Miller, D.L. and Morris, C.D., 2022. Sympatric seals, satellite tracking and protected areas: habitat-based distribution estimates for conservation and management. Frontiers in Marine Science.

Matthiopoulos, J., The use of space by animals as a function of accessibility and preference. Ecological Modelling, 2003. 159(2-3): p. 239-268.

Onoufriou, J., Russell, D.J., Thompson, D., Moss, S.E. and Hastie, G.D., 2021. Quantifying the effects of tidal turbine array operations on the distribution of marine mammals: Implications for collision risk. Renewable Energy, 180, pp.157-165.

Russell, D.J., Hastie, G.D., Thompson, D., Janik, V.M., Hammond, P.S., Scott‐Hayward, L.A., Matthiopoulos, J., Jones, E.L. and McConnell, B.J., 2016. Avoidance of wind farms by harbour seals is limited to pile driving activities. Journal of Applied Ecology, 53(6), pp.1642-1652.

Scott-Hayward, L. 2013. MRSea Package (version 0.1.5): Statistical Modelling of Bird and Cetacean Distributions in Offshore Renewables Development Areas. R Package.

