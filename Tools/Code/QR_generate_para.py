#!/usr/bin/env python3
import sys
import os
import numpy as np


def QR_parafile(gjff, numaH, numa):
    name, tfile=os.path.splitext(gjff)
    fr = open(gjff,"r")
    lines = fr.readlines()
    fr.close()

    ml, sl, b=gjfkeylines(lines)

    #method lines
    mlines=lines[ml:sl[0]]

    #atoms lines
    atomlines=lines[sl[1]+2:sl[2]]
    natoms=len(atomlines)
    
    if natoms != numaH:
        print("Atom number %d in gjf is not equal to %d in pdbH" % (natoms, numaH))
        sys.exit()
    #write parameter file
    fw = open('parameter.dat','w')
    
    fw.write('GAUFILE\n')
    fw.write('%s\n'%gjff)
    fw.write('METHOD\n')
    fw.write('2\n')
    fw.write('IOPT\n')
    fw.write('3\n')
    fw.write('ILINE\n')
    fw.write('1\n')
    fw.write('INITHESSIAN\n')
    fw.write('4\n')
    fw.write('IMULTISTATE\n')
    fw.write('0\n')
    fw.write('UPDATE\n')
    fw.write('3\n')
    
    fw.write('NAT\n')
    fw.write('%d\n'%numa)
    fw.write('HNAT\n')
    fw.write('%d\n'%numaH)
    if os.path.exists('mapfile.dat'):
        fw.write('MAPPING\n')
        fw.write('1\n')
        fw.write('MAPFILE\n')
        fw.write('mapfile.dat\n')
    else:
        fw.write('MAPPING\n')
        fw.write('2\n')
        
    if os.path.exists('res_num.dat'):
        fw.write('RESNUM\n')
        fw.write('1\n')
        fw.write('RESFILE\n')
        fw.write('res_num.dat\n')
    else:
        fw.write('RESNUM\n')
        fw.write('1\n')
    tails = '''IMICROITER
0
MAXENE
10000
MAXSTEP
0.05
LBFGS_MEM
100
ICOORD
0
TOLERANCE
0.0045 0.00001
MAXCYCLE
2000
PRINT
2 2
NCON
0 0
END

'''    
    fw.write(tails)
    fw.close()
    

def QR_datafiles(pdbHf, pdbf):
    pdbatomsH, coordsH, res_numH = Get_pdbinfo(pdbHf)
    pdbatoms, coords, res_num = Get_pdbinfo(pdbf)
    
    maplist = []
    for i, atom in enumerate(pdbatoms):
        try:
            pos = pdbatomsH.index(atom)
            maplist.append(pos+1)


            diff = np.linalg.norm(np.array(coordsH[pos]) - np.array(coords[i]))
            if diff > 0.1:
                print('Warning: Atoms %s in %s far away (%6.3f A) from %s'%(atom,  pdbf, diff, pdbHf))
            
        except ValueError:
            print("Atom %s not found in %s" % (atom, pdbHf))
            sys.exit()
    
    #write mapfile
    fw = open('mapfile.dat','w')
    for i in maplist:
        fw.write("%d\n" % (i))
    fw.close()
    #write res_num file
    fw = open('res_num.dat','w')
    for i in res_numH:
        fw.write("%d\n" % (i))
    fw.close()
    return len(pdbatomsH), len(pdbatoms)

def Get_pdbinfo(pdbf):

    pdbatoms = []
    coords = []
    res_num = []
    chainnum = []
    
    fr = open(pdbf, "r")
    pdbHlines = fr.readlines()
    fr.close()
    

    for line in pdbHlines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            
            keystr = line[11:26]
            value = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            coords.append(value)
            keystr = keystr.replace('WAT', 'HOH')
            keystr = keystr.replace(' OW ', ' O  ')
            keystr = keystr.replace('HID', 'HIS')
            keystr = keystr.replace('HIE', 'HIS')
            keystr = keystr.replace('HIP', 'HIS')
            keystr = keystr.replace('LYN', 'LYS')
            keystr = keystr.replace('GLH', 'GLU')
            keystr = keystr.replace('ASH', 'ASP')
            keystr = keystr.replace('CYM', 'CYS')
            keystr = keystr.replace('CYX', 'CYS')
                
            pdbatoms.append(keystr)
            
            chainstr = line[20:26]
            resnum = int(line[22:26])
            
            if resnum in res_num and chainstr not in chainnum:
                res_num.append(resnum+1000)
            else:
                res_num.append(resnum)
            chainnum.append(chainstr)

    return pdbatoms, coords, res_num

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

if __name__ == '__main__':
    if len(sys.argv) == 4:
        gjfn = sys.argv[1]
        pdbfH = sys.argv[2]
        pdbf = sys.argv[3]
        numaH, numa = QR_datafiles(pdbfH,pdbf)
        QR_parafile(gjfn, numaH, numa)
    else:
        print('Usage: inputs: GauONIOM.gjf PDB_ligH.pdb, mm3.pdb')



