
## Mode of Operation

`aoSystem` has several modes of operation, which determine the calculations made and the format of the output. The mode can be set on the command line with (for example)
```
$ ./aoSystem --mode=C0Var
```
and in the configuration file with
```
mode=C0Var
```
The various choices for mode are described below.

### Error Budgets

If `mode=ErrorBudget` is set, the output is a table with one row per guide star magnitude, showing the optimum actuator pitch, error contributions, and the Strehl.  An example table is:
```
#mag   d_opt  Measurement  Time-delay    Fitting  Chr-Scint-OPD Chr-Index  Disp-Ansio-OPD  NCP-error  Strehl
#mag        d_opt         Measurement      TimeDelay     Fitting         ChrScintOPD        ChrIndex        DispAnisoOPD    NCP              Strehl
0           0.135417      0.0388387        0.0736391     0.277402        0.00646            8.53237e-06     0.0546208       0.287341        0.844109
2.5         0.135417      0.10378          0.0882459     0.277402        0.00646            8.53237e-06     0.0546208       0.287341        0.834352
5           0.135417      0.229569         0.167019      0.277402        0.00646            8.53237e-06     0.0546208       0.287341        0.784162
7           0.270833      0.259455         0.186877      0.50264         0.00081327         8.50977e-06     0.0473598       0.287341        0.644233
8           0.270833      0.352426         0.254413      0.50264         0.00081327         8.50977e-06     0.0473598       0.287341        0.590737
9           0.270833      0.478412         0.346785      0.50264         0.00081327         8.50977e-06     0.0473598       0.287341        0.503288
10          0.40625       0.486257         0.355005      0.697715        0.000239791        8.47918e-06     0.0426013       0.287341        0.393112
10.5        0.40625       0.56633          0.414779      0.697715        0.000239791        8.47918e-06     0.0426013       0.287341        0.345087
11          0.40625       0.659349         0.484968      0.697715        0.000239791        8.47918e-06     0.0426013       0.287341        0.28906
11.5        0.40625       0.767278         0.567574      0.697715        0.000239791        8.47918e-06     0.0426013       0.287341        0.227173
12          0.677083      0.544064         0.409135      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.141631
12.5        0.677083      0.633225         0.478651      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.119886
13          0.677083      0.736571         0.560617      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.095564
13          0.677083      0.736571         0.560617      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.095564
14          0.677083      0.994172         0.772802      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.0461098
14.5        0.677083      1.15307          0.910365      1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.0260051
15          0.677083      1.3354           1.07557       1.18632         3.14861e-05        8.35765e-06     0.0347838       0.287341        0.0118995
```

### The C Tems

These modes calculate the contrast `C` terms due to various contribution as defined in Guyon (2005). The implemented terms are:

- **C0** Variance due to the uncorrected phase, with the effects of scintillation included. 

- **C1** Variance due to the uncorrected amplitude due to scintillation.

- **C2** Variance due to the phase residual (after correction) due to time delay and measurement errors.

- **C3** Variance due to the amplitude residual (after correction) due to time delay and measurement errors.

- **C4** Variance due to the chromaticity of scintillation, causing the ODP measurement to be slightly incorrect at the science wavelength.

- **C5** Variance due to the chromaticity of scintillation, causing the amplitude measurement to be slightly incorrect at the science wavelength.

- **C6** Variance due to the index of refraction of air being wavelength dependent,  causing the ODP measurement to be slightly incorrect at the science wavelength

- **C7** Variance due to atmospheric dispersion causing light at different wavelengths to take different paths through atmospheric turbulence. See Fitzgerald in prep.

 C0 through C6 are as defined in Guyon (2005) as updated in Males & Guyon (2018. C7 is from Fitzgerald in prep. You can choose to calculate a single term, in which case you will get both an ASCII radial profile and 2D FITS image.  An important detail is whether or not a PSF convolution is performed (see Males & Guyon 2018 for details). `aoSystem` will do either, giving the raw variance or the convolved contrast.

NOTE: C3 is not actually implemented, and will return all 0s.  This also means C6 is not implemented and will return all 0s.

| mode      |  Description                                                        | Output
|-----------|---------------------------------------------------------------------|----------------------------|
| C`N`Var   | Calculates the variance contribution from a specific term, where `N=0..7`.  Not convolved.| Radial profile as ASCII 2-columns, separation in [lam/D] and contrast, and FITS containing the contrast in 2D |
| C`N`Con   | Calculates the contrast contribution from a specific term, where `N=0..7`.  This is variance convolved with the PSF.  | Radial profile as ASCII 2-columns, separation in [lam/D] and contrast, and FITS containing the contrast map in 2D|
| CVarAll   | Calculates the variance contributions fom all the terms implemented.  | Radial profiles as ASCII 7-columns, separation in [lam/D] and contrast |
| CConAll  | Calculates the contributions fom all the terms implemented.  This is variance with PSF convolution | Radial profiles as ASCII 7-columns, separation in [lam/D] and contrast |

## Parameters

All parameters are output to a file, which by default is `aoSystem_setup.txt`.  This file can be changed with the `--setupOutFile` CLI option or in the config file with the `setupOutFile` keyword.

### The Atmosphere

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|

### The PSD

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|

### The System

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|
