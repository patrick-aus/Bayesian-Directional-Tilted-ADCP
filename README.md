# Bayesian-Directional-Tilted-ADCP
routines to generate directional spectra from beam-wise velocity measurements from an ADCP

The Bayesian directional method (Hashimoto and Konbune 1988) is used to generate directional spectra from the beam-wise velocity measurements from Nortek Signature1000 ADCP.

The code is presented in MATLAB and builds upon the code presented in Matsuba et al. (2022). The difference here are:
* that no effort is made to reconstruct the infra-gravity waves owing to the non-linear interaction of the sea-swell band waves
* the ***instrument tilt*** is explicity accounted when applying a transfer function



References:
<div class="csl-entry">Hashimoto, N., &#38; Konbune, K. (1988). Directional Spectrum estimation from a Bayesian approach. <i>Coastal Engineering Proceedings</i>, <i>1</i>(21), 4. https://doi.org/10.9753/icce.v21.4</div>
<div class="csl-entry">Matsuba, Y., Roelvink, D., Reniers, A. J. H. M., Rijnsdorp, D. P., &#38; Shimozono, T. (2022). Reconstruction of Directional Spectra of Infragravity Waves. <i>Journal of Geophysical Research: Oceans</i>, <i>127</i>(7), e2021JC018273. https://doi.org/https://doi.org/10.1029/2021JC018273</div>
