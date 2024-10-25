# ONIOM_QR
## Quantum refinement based on ONIOM schemes

### The source code is under construction for more user-friendly.

Quantum refinement which was proposed to combine crystallographic data with computational chemistry methods by several groups can improve the local structures of some proteins. The ONIOM_QR is a quantum refinement method combining several multiscale computational schemes with experimental (X-ray diffraction) information was developed for proteins (especially for metalloproteins and drug-protein systems). 

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

### MLPs
* [ANI](https://aiqm.github.io/torchani/index.html): $\color{red}{Warning:}$ ANI should be run firstly in python cmd to download the parameter files.
* [MLatom](http://mlatom.com/)
* [QDpi](https://gitlab.com/RutgersLBSR/qdpi)

### Tools
Tools and scripts are offered in scripts directory, the source codes are offered in source directory.

## Usage
[dl_gaucns_env_sh](./dl_gaucns_env_sh) should be adjusted with correct environment path including the find.x and QM codes path.

### Workflow
1. PDB preparation
    * Download PDB files (structure & X-ray strcuture factor (.pdb, .cif))
    * Strucutre check (residue flip, missing residue, etc.) e.g: [Saves v6.0](https://saves.mbi.ucla.edu/)
2. QM preparation
    * H add (including protonation state check of residues)
    * ligand/inhibitor H add and FF parameter (RESP: ligand_H.pdb)
    * Protein + ligand/inhibitor (H) check (especially around the ligand/inhibitor): protein_ligH.pdb
    * ONIOM file preparation (GaussianView: protein_ligH.pdb --> protein_ligH.gjf), and H optmizaiton (with heavy atom fixed) using ONIOM(DFT:MM) 
3. CNS preparation
    * [PRODRG](davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg/run.html)/ [ATB](https://atb.uq.edu.au/index.py) using pdb structure
    * using scripts to prepare CNS files (e.g. [dlg_pre](/scripts/dlg_pre),[dlg_cns](/scripts/dlg_cns))
4. QR preparation
    * Put files (protein_ligH.pdb, protein_ligH.gjf and CNS files) in the same directory, using [QR_generate_para.py](/Tools/Code/QR_generate_para.py) to generate mapfile, resfile and parameter.dat files for QR
5. Run QR job
### Examples
Some simple examples are offered including multi-center scheme, quantum refinements using ONIOM scheme and other QM codes.

## Documentation
A simple Manual is offered [Manual_ONIOM_QR.pdf](./Manual_ONIOM_QR.pdf).

## Citations
YAN, Z. Y., Wei, D. C., Li, X., and Chung, L. W. Accelerating Reliable Multiscale Quantum Refinement of Protein–Drug Systems Enabled by Machine Learning. *Nat. Commun.* **2024**, 15, 4181. [DOI: 10.1038/s41467-024-48453-4](https://www.nature.com/articles/s41467-024-48453-4) 

YAN, Z. Y., Li, X., and Chung, L. W. Multiscale Quantum Refinement Approaches for Metalloproteins. *J. Chem. Theory Comput.* **2021**, 17, 3783. [DOI: 10.1021/acs.jctc.1c00148](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00148)

## Special Note to Users
ONIOM_QR is still in the experimental stage and we do not guarantee it will work flawlessly in all your applications.
The micro-iteration in ONIOM is not applied, because of the compatibility between DL-FIND and Gaussian, we are working on this.   



