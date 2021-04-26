# ONIOM_QR
## Quantum refinement based on ONIOM schemes

### The source code is under construction for more user-friendly.

Quantum refinement which was proposed to combine crystallographic data with computational chemistry methods by several groups can improve the local structures of some proteins. The ONIOM_QR is a quantum refinement method combining several multiscale computational schemes with experimental (X-ray diffraction) information was developed for proteins (especially for metalloproteins). 

This code is an interface implemented in an open-source DL-FIND library.
It is developped by [Oscar lab](https://faculty.sustech.edu.cn/oscarchung/en/) in Southern University of Science and Technology .

## Installation
### DL-FIND optimizer
[DL-FIND](https://www.chemshell.org/dl-find) open-source geometry optimisation library should be downloaded.
Then [main_ONIOM_QR.f90](./main_ONIOM_QR.f90) Fortran script is copied to replace the main.f90 file in DL-FIND folder. 
Adjust parameter in the Makefile and complier the code to obtain the executable file find.x.

### CNS
Download and install [CNS](http://cns-online.org/v1.3/), 1.3 version is preferred, other version is not guaranteed.

### QM codes
* [Gaussian](http://gaussian.com/) (ONIOM scheme required)
* [xTB](https://github.com/grimme-lab/xtb) 
* [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal)
* [LSQC](https://itcc.nju.edu.cn/lsqc/)
* [Molpro](https://www.molpro.net/) (WF-in-DFT, high level in ONIOM scheme)
  
### Tools
Tools and scripts are offered in scripts directory, the source codes are offered in source directory.

## Usage
[dl_gaucns_env_sh](./dl_gaucns_env_sh) should be adjusted with correct environment path including the find.x and QM codes path.

### Examples
Some simple examples are offered including multi-center scheme, quantum refinements using ONIOM scheme and other QM codes.


## Citations

