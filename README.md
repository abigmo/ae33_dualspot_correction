# ae33_dualspot_correction
Implementation of the Aethalometer AE33 Dual-spot® correction (Drinovec et al., 2015) applied to MicroAethalometers MA200 (firmware v1.4).
The code is written in R and released alongside a scientific article Bigi et al. (2022) 
It allows to compute the absorption corrected for loading at a custom time step, starting from the raw Attenuation provided by the instrument.

Ideally one can reconstruct also the Attenuation from the instrumental raw count.

Currently the code is applied to 1 minute data for recomputing the corrected absorption at 5 minute intervals.

There are several hard coded parameters which can be changed: the ATN lower and upper limit for FVRF estimate, the Cref and the Mass Absorption Coefficient.

**How to use it**: firstly set your favourite parameters, then load the function `dualspot.compensed.bc` (at the end of the script) and finally compute the data average at your custom time step.

**How to acknowledge this code**:
Only use lower case letters when mentioning `ae33_dualspot_correction`. Always include the version number: ideally, you should also include the Digital Object Identifier (DOI) associated to the specific release you have been using:

ae33_dualspot_correction 

If `ae33_dualspot_correction` was useful for your research, please cite the dedicated article:

Bigi et al., 2022

`ae33_dualspot_correction` relies on external R libraries that require & deserve to be acknowledged in their own right. The following LaTeX blurb is one way to do so:

This research has made use of \textit{ae33_dualspot_correction v1.0.0} \citep[DOI:10.5194/amt-13-6945-2020][]{Bigi_2022}. \textit{ae33_dualspot_correction} relies on the following R packages: \textit{dplyr}, \textit{Rfast}, \textit{tidyr}, \textit{magrittr}


**References**

Drinovec, L., Močnik, G., Zotter, P., Prévôt, A. S. H., Ruckstuhl, C., Coz, E., Rupakheti, M., Sciare, J., Müller, T., Wiedensohler, A., and Hansen, A. D. A.: The "dual-spot" Aethalometer: an improved measurement of aerosol black carbon with real-time loading compensation, Atmos. Meas. Tech., 8, 1965–1979, https://doi.org/10.5194/amt-8-1965-2015, 2015.

