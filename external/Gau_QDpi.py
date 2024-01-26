#!/usr/bin/env python3
import sys
import os
import time
import datetime
import numpy as np

from dpdata_qdpi import QDPiDrivers
from dpdata import System


periodic_table = """ X
  H                                                                                                                           He
  Li  Be                                                                                                  B   C   N   O   F   Ne
  Na  Mg                                                                                                  Al  Si  P   S   Cl  Ar
  K   Ca  Sc                                                          Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr
  Rb  Sr  Y                                                           Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe
  Cs  Ba  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn
  Fr  Ra  Ac  Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr  Rf  Db  Sg  Bh  Hs  Mt  Ds  Rg  Cn  Nh  Fl  Mc  Lv  Ts  Og
""".strip().split()

# QDpi 
'''
长度、能量、距离单位分别是Angstrom，eV，eV/Angstrom，与LAMMPS的metal unit保持一致。

hessian目前不支持。（当然理论上是没有问题的，可能需要讨论一下怎么实现）

'''
def QDpi_run(xyzf, charge, derivatives):
    qdpi = QDPiDriver(
        model="/share/pubbin/qdpi-1.0.pb",
        charge=charge,
        backend="dftb+",
    )
    
    mol = System(xyzf)
    p = mol.predict(driver=qdpi)
    
    ene = p["energies"][0]/27.21138386 #eV 2 bohr
    if derivatives == 0:
        return ene
    elif derivatives == 1:
        forcearray = np.multiply(p["forces"][0], 0.52917720859/27.21138386)
        return ene, forcearray

def mol2xyz(xyzf, eles, coordlist):
    numa = len(eles)
    
    fw = open(xyzf, 'w')
    fw.write('%d\n' % numa)
    fw.write('\n')
    for i in range(numa):
        fw.write('%-16s%14.8f%14.8f%14.8f \n' % (periodic_table[eles[i]], coordlist[i][0], coordlist[i][1], coordlist[i][2]))
    fw.write('\n')
    fw.close()

time_start = time.time()
if len(sys.argv) == 7:
    filein = sys.argv[2]
    fileout = sys.argv[3]
else:
    print('Gaussian external interface using DPpi')
    print('')
    print('DPpi(H C N O elements DFTB+ deltaML)')
    print('')
    print('7 inputs: Gaussian external 6 input')
    print('Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile')
    print('More infomation see Gaussian website')
    sys.exit()

with open(filein, 'r') as f:
    [atoms, deriva, charge ,spin] = [int(s) for s in f.readline().split()]
    ele = np.zeros((atoms,), dtype = np.int)

    atomlist = np.zeros((atoms,3), dtype = np.float32)

    for i in range(atoms):
        axestr = f.readline().split()
        ele[i] = int(axestr[0])
        atomlist[i][0] = float(axestr[1])*0.52917720859
        atomlist[i][1] = float(axestr[2])*0.52917720859
        atomlist[i][2] = float(axestr[3])*0.52917720859

xyzf = os.path.splitext(filein)[0]+'.xyz'
mol2xyz(xyzf, ele, atomlist)

if deriva == 0:
    ene = QDpi_run(xyzf, charge, deriva)
elif deriva == 1:
    ene, force= QDpi_run(xyzf, charge, deriva)
elif deriva == 2:
    #ene, force, hess = QDpi_run(ele, atomlist, deriva)
    sys.exit('Hessian was not supported for now')

time_end = time.time()
timstot=time_end - time_start
print(str(datetime.timedelta(seconds=timstot)))

fw = open(fileout,"w")
fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (ene, 0.0, 0.0, 0.0))
if deriva == 1 : #first derivatives
    for i in range(atoms):       
        fw.write('%20.12E%20.12E%20.12E\n' % (-force[i][0], -force[i][1],-force[i][2]))
fw.close()

