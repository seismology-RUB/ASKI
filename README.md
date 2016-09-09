# ASKI main package

ASKI is a highly modularized program suite offering sensitivity and regularization analysis tools 
for seismic datasets as well as a scattering-integral-type full waveform inversion concept based 
on waveform sensitivity (Fréchet) kernels derived from Born scattering theory (Gauss-Newton 
convergence). ASKI does not implement an intrinsic code for simulation of seismic wave propagation 
but instead comes with support for several external forward codes for 1D and 3D background media 
in spherical and Cartesian framework, at the moment 
[SPECFEM3D_Cartesian](https://github.com/seismology-RUB/SPECFEM3D_Cartesian_for_ASKI), 
[SPECFEM3D_GLOBE](https://github.com/seismology-RUB/SPECFEM3D_GLOBE_for_ASKI), 
[Gemini II](https://www.geophysik.rub.de/trac/gemini) and 
[NEXD](http://www.rub.de/nexd). 


## Authors and License

ASKI and some of its components, as well as documentation and some examples
are available under terms of the [GNU General Public License](LICENSE) (version 2 or higher)
on [github](https://github.com/seismology-RUB/ASKI).
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/aski in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the [github repository](https://github.com/seismology-RUB/ASKI).

The main authors are Florian Schumacher and Wolfgang Friederich (Ruhr-University Bochum, Germany).


## Documentation

Please refer to documents in [doc/](doc/) :

* the [ASKI user manual](doc/ASKI_manual.pdf)
* [Florian Schumacher's doctoral dissteration](doc/dissertation_florian_schumacher.pdf) 
  about waveform sensitivity kernels and the modularized iterative full waveform inversion 
  concept on which ASKI is based
* the accepted version of our [GJI 2016 paper](doc/ASKI_paper_gji_2016.pdf), re-typeset in 
  a standard layout


## Toy examples

Any files packages for the ASKI toy examples, as described in the [manual](doc/ASKI_manual.pdf) (chapter 0 *ASKI workflows*) are attached to [release 1.0 of the ASKI main package](https://github.com/seismology-RUB/ASKI/releases/tag/v1.0). These inversion examples use release version 1.0 of the extension packages [SPECFEM3D_Cartesian](https://github.com/seismology-RUB/SPECFEM3D_Cartesian_for_ASKI/releases/tag/v1.0) and [SPECFEM3D_GLOBE](https://github.com/seismology-RUB/SPECFEM3D_GLOBE_for_ASKI/releases/tag/v1.0), respectively.


## Requirements

* GNU Make
* Fortran compiler (sufficient standard, so far tested: GNU Fortran 
  (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609)
* BLAS and LAPACK  libraries for all applications
* BLACKS, SCALAPACK and MPI libraries optionally for few parallel applications 
  (e.g. available via [netlib.org/](http://www.netlib.org/))


## Installation

0. If not yet done, you should downloade the source code of the ASKI main package by either
   * cloning the master branch of the ASKI repository on gitHub.com:
     ```
     git clone --depth 1 --branch master https://github.com/seismology-RUB/ASKI
     ```
     
   * or downloading a zipped version of the source code from [there](https://github.com/seismology-RUB/ASKI/archive/master.zip):
     ```
     wget https://github.com/seismology-RUB/ASKI/archive/master.zip
     ```
     
     The directory to where you have cloned/extracted the master branch, will in the following be referred to as `ASKI/`
1. Adjust the software to your system and personal requirements by:
   * setting the variables `COMPILER`, `MPICOMPILER` in file [ASKI/Makefile](Makefile) appropriately
   * setting the variables `BLAS`, `LAPACK`, (`BLACS`, `SCALAPACK`, `MPILIB`) in file [ASKI/Makefile](Makefile) in 
     order to correctly bind required libraries
   * adjusting the variable `FFLAGS` in file [ASKI/Makefile](Makefile), if required

2. Run command
   ```
   make all
   ```
   from installation path `ASKI/` to compile all serial programs (nearly all functionality)

3. Run command
   ```
   make parallel
   ```
   from installation path `ASKI/` to compile all parallel programs (two optional linear system solvers which run in parallel)

Now [ASKI/bin](bin/) should contain all binaries and no error should have occurred (at best...).


## Usage

Refer to the [user manual](doc/ASKI_manual.pdf) for any details on using ASKI.
