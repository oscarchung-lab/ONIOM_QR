# Tools for QR works
## 

### Amber2GauFF.py
    generate formce field in Gaussian format (.prm) from Amber force field
    python3 Amber2GauFF.py amber.dat

### FFcombine.sh
    Combine two Amber force parameter files (e.g. parm10.dat + frcmod.ff14SB)

### PDB2pqr.sh
    call PDB2PQR(https://server.poissonboltzmann.org/pdb2pqr_2.0.0/) 

### Check_h2o_Connect.sh
    check connectivity of H2O, sometimes two O are too closed, the GV gives wrong connectivity.

### GauOniomchr.sh
    Get total charge of H, M, L layer

### GauOniomLink.sh
    Set relaxed atoms (where the link atoms should be relaxed)

### Check_h2o_Connect.sh
    Reoder the connectivity of water in the gjf

### RESP_mol2gjf.py
    translate the RESP file (.mol2 file from anterchamber) to gjf format (high layer with relaxed all toms in ONIOM)

### AB_Conformer_V1.f90
    separate comformer of PDB to A and B files

### QR_High.f90
    This program is to modifiy the layer atom of gjf file
    Input: mol.gjf, High.xyz H (H or M)'

### QGetRSZD_res.sh
    This program is to extract the RSZD scores from EDSTATS log.txt
    Input: PDB

