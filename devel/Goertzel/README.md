## executable `xgoertzel`

Program [goertzel.f90](goertzel.f90) is an example of how to implement discrete Fourier
transforms of a signal in an *on-the-fly* fashion, as might be required
by time-stepping forward codes that are to be extended for ASKI.

The program reads in the synthetic seismogram file 
[data_SB04_OUTPUT_FILES--NETW.RC02.FXY.semd](data_SB04_OUTPUT_FILES--NETW.RC02.FXY.semd)
that was used as "measured data" of event SB04, station RC02, component CY
in the "cross borehole" example inversion provided on 
[gitHub](https://github.com/seismology-RUB/ASKI/releases/tag/v1.0)
(in particular, file 
`ASKI_inversion_cross_borehole/measured_data/data_SB04_OUTPUT_FILES/NETW.RC02.FXY.semd`).
It is also shown in the [GJI 2016 ASKI paper](http://dx.doi.org/10.1093/gji/ggv505) 
Figure 7, upper-right, and in 
[Florian's dissertation](http://nbn-resolving.de/urn/resolver.pl?urn=urn:nbn:de:hbz:294-40511) 
Figure 5.4 (b).

With the same time series, the program performs two different implementations of
discrete Fourier transform to the discrete set of frequencies as used in the
"cross borehole" example inversion:
1. "conventional" explicit summation
2. Goertzel's algorithm in "reversed" form, i.e.\ correcting for time-reversal and phase shift
   (see [ASKI developers manual](../doc/ASKI_developers_manual.pdf), section 4.5)

Moreover, `xgoertzel` measures the performance of the time loop in both algorithms (by repeating
it for 50000 times, please adapt parameter variable `number_of_repetitions` in [goertzel.f90](goertzel.f90) 
as required). Other computational resources (pre-processing such as pre-computing `efactors`, 
memory requirements, as well as post-processing such as phase-shifting Goertzel's result) is *not* 
accounted for in the performance test. Using compiler flags `-O1`, `-O2` or `-O3` (and *not* enabling 
debug compiler flags such as `-fbounds-check`, `-fbacktrace`), Goertzel's algorithm can be more 
than 2 times faster than the conventional explicit summation.


### compile and run

If required, edit [Makefile](Makefile), in particular values of `COMPILER` and `FFLAGS`.

Compile and immediately run the executable by
```
make
```

or 
```
make all
```

Compile the executable `xgoertzel` only (and not run it) by
```
make compile
```

Run the compiled executable `xgoertzel` only (without compiling it beforehand) by
```
make run
```

or simply
```
./xgoertzel
```

### program output

The program's output is on screen only, no output files are written. For plotting the output
spectra, the screen output can be copied to text files. This was already done here for
[conventional summation](out_conventional-DFT.txt), 
[original Goertzel's algorithm](out_Goertzel-DFT-original.txt) (code section is commented
for benchmark) and for [adapted Goertzel's algorithm](out_Goertzel-DFT-reverse.txt).

For comparison, the respective transformed [measured data file](measured_data--data_SB04_RC02_CY)
of the "cross borehole" example was copied here from 
`ASKI_inversion_cross_borehole/measured_data/data_SB04_RC02_CY`. Please note, that this file was 
produced by executable `transformSpecfem3dCartesianMeasuredData` of the 
[SPECFEM3D for ASKI](https://github.com/seismology-RUB/SPECFEM3D_Cartesian_for_ASKI) package.
This executable uses module [discreteFourierTransform](../../f90/discreteFourierTransform.f90), 
which in principle implements an explicit summation in the same way as done by program 
[goertzel.f90](goertzel.f90). The slight numerical differences in the output spectra might stem from 
an older compiler version / older ASKI version used to produce file `measured_data--data_SB04_RC02_CY`.

All four spectrum files can be plotted running the script [plot_spectra.sh](plot_spectra.sh).
