#!/usr/bin/env python3

import sys
import os
import numpy as np
import re

eleslist = ['H','He',
        'Li','Be','B','C','N','O','F','Ne',
        'Na','Mg','Al','Si','P','S','Cl','Ar',
        'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']

def run_Mlatom(method, filein, fileout):
    ele, atomlist, deriva, charge, spin = get_external_coord(filein)

    if deriva == 0: #ene
        ene = Mlatom_run(method, ele, atomlist, deriva, charge, spin)
        fw = open(fileout,"w")
        fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (ene, 0.0, 0.0, 0.0))
        for i in range(len(ele)):       
            fw.write('%20.12E%20.12E%20.12E\n' % (0.0, 0.0, 0.0))
        fw.close()
    elif deriva == 1: #ene + grad 
        ene, grad = Mlatom_run(method, ele, atomlist, deriva, charge, spin)
        fw = open(fileout,"w")
        fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (ene, 0.0, 0.0, 0.0))
        for i in range(len(ele)):       
            fw.write('%20.12E%20.12E%20.12E\n' % (grad[i][0], grad[i][1],grad[i][2]))
        fw.close()

def get_external_coord(filein):
    with open(filein, 'r') as f:
        [atoms, deriva, charge ,spin] = [int(s) for s in f.readline().split()]
        ele = np.zeros((atoms,), dtype = int)
        
        atomlist = np.zeros((atoms,3), dtype = float)

        for i in range(atoms):
            axestr = f.readline().split()
            ele[i] = int(axestr[0])
            atomlist[i][0] = float(axestr[1])*0.52917724
            atomlist[i][1] = float(axestr[2])*0.52917724
            atomlist[i][2] = float(axestr[3])*0.52917724
    return ele, atomlist, deriva, charge, spin

'''
     imult   Definition of multiplicity.
              *** Options for RHF calculations.
              = 0 Closed-shell singlet.
              = 1 Open-shell singlet with two singly occupied orbitals.
                  This usually corresponds to an excited singlet state.
              = 2 Doublet.
              = 3 Triplet.
              *** Options for UHF calculations.
              = 0 Singlet.
              = 1 Singlet (same as imult=0).
              = 2 Doublet.
              = 3 Triplet.
              = 4 Quartet (etc).
            *** Note on RHF and UHF calculations.
              By default, RHF for imult=0 and UHF for imult.gt.0
              which may be changed using option iuhf (see below).
              imult.gt.3 is possible for UHF only.
'''
def Mlatom_run(method, eles, coordlist, derivatives, chg, spin):
    if chg != 0 or spin !=1 :
        print('charge = %5d and mult = %5d' %(chg, spin))
        if not os.path.exists('mndokw'):     
            fw = open('mndokw',"w")
            fw.write('iop=-22 immdp=-1 + \n')
            fw.write('igeom=1 iform=1 + \n')
            fw.write('jop=2 nsav15=3 + \n')
            fw.write('iprint=-1 kprint=-5 lprint=-2 mprint=0 jprint=-1 + \n')
            fw.write('kharge=%-3d imult=0 nprint=-1  \n'%(chg))
            fw.close()
            print('A mndokw file have been created!')
    if not os.path.exists('gau_mlatom.inp'):
        fw = open('gau_mlatom.inp',"w")
        fw.write(method+' \n')
        fw.write('xyzfile=gau_mlatom.xyz \n')
        fw.write('yestfile=ene_mlatom.dat \n')
        if derivatives == 1 :
            fw.write('ygradxyzestfile=grad_mlatom.dat \n')
        if chg != 0 or spin !=1 :
            fw.write('mndokeywords=mndokw \n')
        fw.close()

    #write xyz file
    fw = open('gau_mlatom.xyz',"w")
    fw.write(str(len(eles)).replace(' ','')+'\n')
    fw.write('gau_mlatom\n')
    for i in range(len(eles)):
        fw.write('%-2s %20.12f%20.12f%20.12f \n'%(eleslist[eles[i]-1], coordlist[i][0], coordlist[i][1], coordlist[i][2]))
    fw.close()
    #run job
    print('Mlatom lanch ...')
    os.system('mlatom gau_mlatom.inp > gau_mlatom.out')
    print('Mlatom job finished!')
    #read energy
    fr = open('ene_mlatom.dat',"r")
    lineene = fr.readline()
    ene = float(lineene)
    fr.close()

    #read grad
    if derivatives == 1 :
        grad = np.zeros((len(eles), 3), dtype = np.float64)
        fr = open('grad_mlatom.dat',"r")
        linetmp = fr.readline()
        linetmp = fr.readline()
        for i in range(len(eles)):
            linetmp = fr.readline()
            vartmp = linetmp.split()
            grad[i][0] = float(vartmp[0])*0.52917721092
            grad[i][1] = float(vartmp[1])*0.52917721092
            grad[i][2] = float(vartmp[2])*0.52917721092
        fr.close()

    os.system('rm -f mndo.* fort.4 fort.11 fort.15 std')
    os.system('rm -f gau_mlatom.xyz ene_mlatom.dat grad_mlatom.dat gau_mlatom.out')
    
    if derivatives == 0:
        return ene
    elif derivatives == 1:
        return ene, grad

if __name__ == '__main__':
    if len(sys.argv) == 7:
        filein = sys.argv[2]
        fileout = sys.argv[3]
        method = 'AIQM1'
    elif len(sys.argv) == 8:
        method = sys.argv[1]
        filein = sys.argv[3]
        fileout = sys.argv[4]
    else:
        print('Gaussian external interface using Mlatom')
        print('Option 1: no input : AIQM1 without mndokw' )
        print('Option 2: AIQM1 : AIQM1 without mndokw')
        print('AIQM1 --> AIQM1@DFT or AIQM1@DFT* please see Mlatom manual')
        print('mndokw: --> mndokeywords file name including charge and multi info')
        print('if charge and mult are not 0, then a simple mndokw will be generated')
        print('Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile')
        print('More infomation see Gaussian website')
        sys.exit()

    run_Mlatom(method, filein, fileout)



