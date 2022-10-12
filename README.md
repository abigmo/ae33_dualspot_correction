# ae33_dualspot_correction
Implementation of the Aethalometer AE33 Dual-spot® correction (Drinovec et al., 2015) applied to MicroAethalometers MA200 (firmware v1.4).
The code is in R and it is provided as is.
It allows to compute the absorption corrected for loading at a custom time step, starting from the raw Attenuation provided by the instrument.

Ideally one can reconstruct also the Attenuation from the instrumental raw count.

Currently the code is applied to 1 minute data and recomputes the corrected absorption at 5 minute intervals.

There are several hard coded parameters which can be changed: the ATN lower and upper limit for FVRF estimate, the Cref and the Mass Absorption Coefficient.

How to use it: firstly set your favourite parameters, then load the function dualspot.compensed.bc (at the end of the script) and finally compute the data average at your custom time step.
