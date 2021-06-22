
## Mode of Operation

`aoSystem` has several modes of operation, which determine the calculations made and the format of the output. The mode can be set on the command line with (for example)
```
$ ./aoSystem --mode=C0Raw
```
and in the configuration file with
```
mode=C0Raw
```
The various choices for mode are described below.

### Error Budgets

If `mode=ErrorBudget` is set, the output is a table with one row per guide star magnitude, showing the optimum actuator pitch, error contributions, and the Strehl.  An example table is:
```
#mag   d_opt  Measurement  Time-delay    Fitting  Chr-Scint-OPD Chr-Index  Disp-Ansio-OPD  NCP-error  Strehl
0      0.135   1.60746      16.2009      23.3996       0            0           0            0        0.951111
2      0.135   4.03782      16.2009      23.3996       0            0           0            0        0.950306
4      0.135   9.59609      16.4651      23.3996       0            0           0            0        0.945371
6      0.135   18.1457      19.6898      23.3996       0            0           0            0        0.924964
8      0.135   32.6527      28.6758      23.3996       0            0           0            0        0.860478
10     0.135   58.905       47.1755      23.3996       0            0           0            0        0.680388
12     0.135   107.522      82.6173      23.3996       0            0           0            0        0.311002
13     0.135   146.079      111.393      23.3996       0            0           0            0        0.120574
14     0.135   199.615      152.385      23.3996       0            0           0            0        0.0197602
14.5   0.135   234.129      179.432      23.3996       0            0           0            0        0.00451157
15     0.135   275.513      212.431      23.3996       0            0           0            0        0.00055321
15.5   0.135   325.649      253.079      23.3996       0            0           0            0        2.68251e-05
16     0.135   387.173      303.664      23.3996       0            0           0            0        3.15645e-07
```

### The C Tems

These modes calculate the contrast of each term C0, C1, C2, C4, C6, and C7. C0 through C4 as defined in Guyon (2005) as updated in Males & Guyon (2018). C6 and C7 are the dispersive anisoplanatism terms from Fitzgerald in prep.  Note that C3 and C5 are not implemented.

You can choose to calculate a single term, in which case you will get both an ASCII radial profile and 2D FITS image.  An important detail is whether or not a PSF convolution is performed (see Males & Guyon 2018 for details). `aoSystem` will do either. 

| mode      |  Description                                                        | Output
|-----------|---------------------------------------------------------------------|----------------------------|
| C<N>Raw   | Calculates the contrast contribution from a specific term, where N=0,1,2,4,6 or 7.  Not convolved.| ASCII 2-column, separation in [lam/D] and contrast, and FITS containing the contrast in 2D |
| C<N>Map   | Calculates the contrast contribution from a specific term, where N=0,1,2,4,6 or 7.  Convolved with the PSF.  | ASCII 2-column, separation in [lam/D] and contrast, and FITS containing the contrast map in 2D|
| CAllRaw   | Calculates the contributions fom all the terms implemented.  | ASCII 7-columns, separation in [lam/D] and contrast |
| CProfAll  | Calculates the contributions fom all the terms implemented, with PSF convolution | ASCII 7-columns, separation in [lam/D] and contrast |

## Parameters

All parameters are output to a file, which by default is `mxAOAnalysisSetup.txt`.  This file can be changed with the `--setupOutFile` CLI option or in the config file with the `setupOutFile` keyword.

### The Atmosphere

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|

### The PSD

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|

### The System

| Param   | CLI      | section.key | units | notes      |
|---------|----------|-------------|-------|------------------|
