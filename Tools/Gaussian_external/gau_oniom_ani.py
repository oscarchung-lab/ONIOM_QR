#!/usr/bin/env python3

import sys
import os
import numpy as np
import torch
import torchani
import re

eleslist = ['H','He',
        'Li','Be','B','C','N','O','F','Ne',
        'Na','Mg','Al','Si','P','S','Cl','Ar',
        'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']
scalematrix=np.loadtxt('/share/home/chem-zhonglh/yanzy/bin/ONIOM_scale.txt',delimiter=',',dtype = np.float64)


def run_Oniom_ANI(gjff, hm, lm, filein, fileout):
    name,tfile=os.path.splitext(gjff)
    fr = open(gjff,"r")
    lines = fr.readlines()
    fr.close()
    ml, sl, ifconn = gjfkeylines(lines)

    #method lines
    mlines = lines[ml:sl[0]]

    #Lchg Lspin Hchg Hspin Hchg Hspin
    chgspin = [int(s) for s in lines[sl[1]+1].split()]
    #atoms lines
    atomlines = lines[sl[1]+2:sl[2]]
    natoms = len(atomlines)
    eles, coords, Hlist, llist = getcoords_frag(atomlines)

    elesnum = np.zeros(natoms, dtype=int)
    for i in range(natoms):
        elesnum[i] = ele2num(eles[i])

    elein, coordsin, deriva = get_external_coord(filein)

    if (elesnum != elein).all():
        sys.exit('The element list in '+gjff+' is different with '+filein)

    if ifconn:
        #generate connectivity matrix from connectivity block
        if len(sl) == 3:
            connlines = lines[sl[2]+1:]
        else:
            connlines = lines[sl[2]+1:sl[3]]

        matrix_link = linkmatrix(connlines)
        linklist = Getlayerlist(matrix_link, Hlist)
    else:
        sys.exit('Please give geom=connectivity and connect info in the '+gjff)
        

    Hdim = len(Hlist)

    Hlistele, Hlistcoord = get_hele_hcoords(eles,coordsin,Hlist,linklist)
    
    # run computations
    if deriva == 0:
        if hm == 'mlatom':
            HM_E = Mlatom_run(Hlistele, Hlistcoord, deriva, chgspin[2], chgspin[3])
            print('High level of Model using mlatom '+hm)
            print('E_HM('+hm+') =', HM_E, '(bohr)')
        else:
            HM_E = ANI_run(Hlistele, Hlistcoord, hm, deriva)
            print('High level of Model using ANI '+hm)
            print('E_HM('+hm+') =', HM_E, '(bohr)')

        if lm == 'mlatom':
            LM_E = Mlatom_run(Hlistele, Hlistcoord, lm, deriva, chgspin[2], chgspin[3])
            print('Low level of Model using mlatom '+lm)
            print('E_LM('+hm+') =', LM_E, '(bohr)')
            LS_E = Mlatom_run(elein, coordsin, lm, deriva, chgspin[0], chgspin[1])
            print('Low level of system using mlatom '+lm)
            print('E_LS('+hm+') =', LS_E, '(bohr)')
        else:
            LM_E = ANI_run(Hlistele, Hlistcoord, lm, deriva)
            print('Low level of Model using ANI '+lm)
            print('E_LM('+hm+') =', LM_E, '(bohr)')
            LS_E = ANI_run(elein, coordsin, lm, deriva)
            print('Low level of system using ANI '+lm)
            print('E_LS('+hm+') =', LS_E, '(bohr)')

        #
        E_ANI = LS_E + HM_E - LM_E
        fw = open(fileout,"w")
        fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (E_ANI, 0.0, 0.0, 0.0))
        fw.close()   
    elif deriva == 1:
        if hm == 'mlatom':
            HM_E, HM_F = Mlatom_run(Hlistele, Hlistcoord, deriva, chgspin[2], chgspin[3])
            print('High level of Model using mlatom '+hm)
            print('E_HM('+hm+') =', HM_E, '(bohr)')
        else:
            HM_E, HM_F = ANI_run(Hlistele, Hlistcoord, hm, deriva)
            print('High level of Model using ANI '+hm)
            print('E_HM('+hm+') =', HM_E, '(bohr)')
            
        if lm == 'mlatom':
            LM_E, LM_F = Mlatom_run(Hlistele, Hlistcoord, lm, deriva, chgspin[2], chgspin[3])
            print('Low level of Model using mlatom '+lm)
            print('E_LM('+hm+') =', LM_E, '(bohr)')
            LS_E, LS_F = Mlatom_run(elein, coordsin, lm, deriva, chgspin[0], chgspin[1])
            print('Low level of system using mlatom '+lm)
            print('E_LS('+hm+') =', LS_E, '(bohr)')
        else:
            LM_E, LM_F = ANI_run(Hlistele, Hlistcoord, lm, deriva)
            print('Low level of Model using ANI '+lm)
            print('E_LM('+hm+') =', LM_E, '(bohr)')
            LS_E, LS_F = ANI_run(elein, coordsin, lm, deriva)
            print('Low level of system using ANI '+lm)
            print('E_LS('+hm+') =', LS_E, '(bohr)')

        F_ANI = LS_F
        #
        E_ANI = LS_E + HM_E - LM_E
        ldimtmp = 0 
        for i in range(Hdim):
            F_ANI[Hlist[i]][:] = LS_F[Hlist[i]][:] - LM_F[i][:] + HM_F[i][:]
            if len(linklist[i]) != 0: #pure high level
                for j, linktmp in enumerate(linklist[i]):
                    scaletmp = scalematrix[eleslist.index(eles[linktmp])][eleslist.index(eles[Hlist[i]])]
                    F_ANI[Hlist[i]][:] = F_ANI[Hlist[i]][:] + (1-scaletmp)*(HM_F[Hdim+j][:]- LM_F[Hdim+j][:])
                    F_ANI[linktmp][:] = F_ANI[linktmp][:] + scaletmp*(HM_F[Hdim+j][:]- LM_F[Hdim+j][:])
                ldimtmp = ldimtmp + len(linklist[i])

        fw = open(fileout,"w")
        fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (E_ANI, 0.0, 0.0, 0.0))
        for i in range(natoms):       
            fw.write('%20.12E%20.12E%20.12E\n' % (F_ANI[i][0], F_ANI[i][1],F_ANI[i][2]))
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
    return ele, atomlist, deriva

def get_hele_hcoords(eles,coords,Hlist,linklist): 
    #Get Model: eles & coords
    Hdim = len(Hlist)
    hele = np.zeros(Hdim, dtype=int)
    hcoord = np.zeros((Hdim, 3), dtype = np.float64)
    lcoord = np.zeros((Hdim, 3), dtype = np.float64)
    lele = np.zeros(Hdim, dtype=int)
    ldim = 0
    for i in range(Hdim):
        hele[i] = ele2num(eles[Hlist[i]])
        hcoord[i][:] = coords[Hlist[i]][:]
        if len(linklist[i]) != 0:
            for j, linktmp in enumerate(linklist[i]):
                lele[ldim+j] = 1
                #lcoord[ldim+j][:] = coords[linktmp][:]
                scaletmp = scalematrix[eleslist.index(eles[linktmp])][eleslist.index(eles[Hlist[i]])]
                lcoord[ldim+j][:] = link2H(coords[Hlist[i]][:], coords[linktmp][:], scaletmp)
            ldim = ldim + len(linklist[i])
    print('There are %5d atoms in high layer, and %5d link atoms!\n' %(Hdim, ldim))

    Modeldim = Hdim+ldim
    Hlistele = np.zeros(Modeldim, dtype=int)
    Hlistcoord = np.zeros((Modeldim, 3), dtype = np.float64)
    Hlistele[:Hdim] = hele[:]
    Hlistele[Hdim:] = lele[:ldim]

    Hlistcoord[:Hdim][:] = hcoord
    Hlistcoord[Hdim:][:] = lcoord[:ldim][:]

    return Hlistele, Hlistcoord


def ele2num(ele):
    return eleslist.index(ele)+1

def Getlayerlist(linkmatrix, fraglist):
    # each frag contain defined by a fraglist
    linkatomlist = []
 
    for atomtmp in fraglist:   
        linktmp = np.nonzero(linkmatrix[atomtmp][:])
        link_atomj = []
        for linkatomtmp in linktmp[0]:
            if linkatomtmp not in fraglist:
                link_atomj.append(linkatomtmp)
        linkatomlist.append(link_atomj)

    return linkatomlist

#get key lines from Gaussian gjf file
def gjfkeylines(lines):
    spacelist=[]
    for i in range(len(lines)):
        #method lines
        if lines[i].startswith('#'):
            mline=i
        #empty lines
        if lines[i].isspace() :
            #repeat empty lines at the end of files
            if len(spacelist)> 1 and i==spacelist[-1]+1:
                break
            spacelist.append(i) 
    #if contains connectivity key word
    ifconn=False
    for linestr in lines[mline:spacelist[0]]:
        if 'geom=connectivity' in linestr:
            ifconn=True
    return mline, spacelist, ifconn

#get coords from atom block
def getcoords_frag(lines):
    natoms=len(lines)
    coords=np.zeros((natoms,3),dtype = np.float64)
    elelist=[]
    Hlist = []
    linkpaire=[]
    linklist = []
    # ele x y z
    for i, linestr in enumerate(lines):
        linktmp = []
        linklist.append(linktmp)
        vartmp=linestr.split()
        elelist.append(vartmp[0])
        coords[i][0]=float(vartmp[1])
        coords[i][1]=float(vartmp[2])
        coords[i][2]=float(vartmp[3])
        if len(vartmp) == 5 and vartmp[4] == 'H':
            Hlist.append(i)
        elif len(vartmp) == 7 and vartmp[4] == 'L':
            linktmp.append(i)
            linktmp.append(int(vartmp[6])-1)
            linkpaire.append(linktmp)
        elif len(vartmp) == 5 and vartmp[4] == 'L':
            continue
        else:
            sys.exit(linestr+' This line is wrong')
    #treat link
    for i, sublist in enumerate(linkpaire):
        if len(sublist) > 0:
            linklist[Hlist.index(sublist[1])].append(sublist[0])
    return elelist, coords, Hlist, linklist

# dis(coord1--H)=scale*(coord1--coord2)
# H replace coord2
def link2H(coord1,coord2,scale):
    return coord1 + scale*(coord2-coord1)

#get link matrix from connectivity block
def linkmatrix(lines):
    linkm=np.zeros((len(lines),len(lines)), dtype = np.float64)
    for i, linestr in enumerate(lines):
        var=linestr.split()
        if len(var) == 1:
            continue
        else:
            j=1
            while j < len(var):        
                linkm[i][int(var[j])-1]=float(var[j+1])
                linkm[int(var[j])-1][i]=float(var[j+1])
                j=j+2
    return linkm

#ANI run
def ANI_run(eles, coordlist, method, derivatives):
    # ANI
    #device = torch.device('cpu')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if method == '1x':
        #print('Method '+method+' was used')
        model = torchani.models.ANI1x(periodic_table_index=True).to(device).double()
    elif method == '1ccx':
        #print('Method '+method+' was used')
        model = torchani.models.ANI1ccx(periodic_table_index=True).to(device).double()
    elif method == '2x':
        #print('Method '+method+' was used')
        model = torchani.models.ANI2x(periodic_table_index=True).to(device).double()
    else: 
        print('Method '+method+' was unkown')
        sys.exit()
    coordinates = torch.from_numpy(coordlist).requires_grad_(True).unsqueeze(0)
    species = torch.from_numpy(eles).unsqueeze(0)
    masses = torchani.utils.get_atomic_masses(species)

    # Now let's compute energy and force:
    energy = model((species, coordinates)).energies
    if derivatives == 2: #first derivatives
        hessian = torchani.utils.hessian(coordinates, energies=energy)
        hess_mat = hessian.numpy()[0] * 0.52917721092 * 0.52917721092
    
    if derivatives == 1 or derivatives == 2: #first derivatives
        derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
        force = derivative
        forcearray=force.squeeze().numpy() * 0.52917721092


    #freq, modes, fconstants, rmasses = torchani.utils.vibrational_analysis(masses, hessian, mode_type='MDU')
    #torch.set_printoptions(precision=12, sci_mode=False)
    
    if derivatives == 0:
        return energy.item()
    elif derivatives == 1:
        return energy.item(), forcearray
    elif derivatives == 2:
        return energy.item(), forcearray, hess_mat

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

def Mlatom_run(eles, coordlist, derivatives, chg, spin):
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
        fw.write('AIQM1 \n')
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
    if len(sys.argv) == 8:
        gjff = sys.argv[1]
        filein = sys.argv[3]
        fileout = sys.argv[4]
        hm='1ccx'
        lm='2x'
    elif len(sys.argv) == 10:
        gjff = sys.argv[1]
        hm = sys.argv[2]
        lm = sys.argv[3]
        filein = sys.argv[5]
        fileout = sys.argv[6]
    else:
        print('Gaussian external interface using ONIOM-like ANI embbeding')
        print('including: Mlatom ANI-1x, ANI-1ccx, ANI-2x--1x, 1ccx, 2x')
        print('ANI-1x(H C N O elements wB97X/6-31G(d))')
        print('ANI-1ccx(H C N O elements CCSD(T)*/CBS (CCSD(T))')
        print('ANI-2x(H C N O F S Cl elements wB97X/6-31G(d)) default')
        print('')
        print('1 input: model.gjf; no extral input: 1ccx:2x ')
        print('3 extral inputs: model.gjf 1ccx 2x <==> 1ccx:2x or mlatom 2x <==> mlatom:2x')
        print('Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile')
        print('More infomation see Gaussian website')
        sys.exit()

    run_Oniom_ANI(gjff, hm, lm, filein, fileout)

