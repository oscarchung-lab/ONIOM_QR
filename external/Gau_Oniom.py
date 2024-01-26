#!/usr/bin/env python3
import sys
import os
import numpy as np

def GaussianOniom(gjff, filein, fileout):
    name, tfile=os.path.splitext(gjff)
    elein, coordsin, deriva = get_external_coord(filein)
    
    numatom = len(elein)
    #update gjff 
    UpdateOniomcoord(gjff, coordsin)
    
    #Run gaussian 
    os.system("g16 < "+gjff+" > "+name+'.log')
    #read energy & gradient
    energy = GetGaussianlogenergy(name+'.log')
    grad = Getfort7frad('fort.7', numatom)

    #write outfile 
    fw = open(fileout,"w")
    fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (energy, 0.0, 0.0, 0.0))
    if deriva == 1 or deriva == 2: #first derivatives
        for i in range(numatom):       
            fw.write('%20.12E%20.12E%20.12E\n' % (grad[i][0], grad[i][1],grad[i][2]))
        #if deriva == 2: #first derivatives
        #    fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))
        #    fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))
        #    for i in range(3*atoms):    
        #        fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))
#
        #    hess_lt = hess[np.tril_indices(3*atoms)]
        #    k = int(len(hess_lt)/3)
        #    for i in range(k):
        #        a = hess_lt[i*3]
        #        b = hess_lt[i*3+1]
        #        c = hess_lt[i*3+2]
        #        fw.write('%20.12E%20.12E%20.12E\n' % (hess_lt[i*3], hess_lt[i*3+1],hess_lt[i*3+2]))
    fw.close()

def Getfort7frad(fort7f, numa):
    fr = open(fort7f,"r")
    lines = fr.readlines()
    fr.close()
    # line 0-numa coordinates
    grad = np.zeros((numa,3), dtype = float)
    # line numa - numa*2
    for i, line in enumerate(lines[numa:numa*2]):
        axestr = line.split()
        grad[i][0] = float(axestr[0].replace("D","E"))
        grad[i][1] = float(axestr[1].replace("D","E"))
        grad[i][2] = float(axestr[2].replace("D","E"))
    return grad
    
def GetGaussianlogenergy(logf):
    fr = open(logf,"r")
    lines = fr.readlines()
    fr.close()
    
    SCFenergy = 0.0 
    Oniomenergy = 0.0
    for line in lines:
        if "SCF Done:" in line:
            SCFenergy=float(line.split()[4])
        if 'extrapolated energy' in line:
            Oniomenergy = float(line.split()[4])
    if Oniomenergy == 0.0:
        return SCFenergy
    else:
        return Oniomenergy
    
    
def UpdateOniomcoord(fgjf, coords):
    name, tfile=os.path.splitext(fgjf)
    fr = open(fgjf,"r")
    lines = fr.readlines()
    fr.close()

    ml, sl, b=gjfkeylines(lines)

    #method lines
    mlines=lines[ml:sl[0]]

    #atoms lines
    atomlines=lines[sl[1]+2:sl[2]]
    natoms=len(atomlines)
    
    #os.rename(fgjf, name+'_ori'+tfile)
    
    fw = open(fgjf,"w")
    fw.writelines(lines[:sl[1]+2])
    #atomlines
    for i, line in enumerate(atomlines):
        xyz = coords[i][:3]
        fw.write('%19s%14.8f%14.8f%14.8f%s'%(line[:19],xyz[0],xyz[1],xyz[2],line[61:]))
    
    fw.writelines(lines[sl[2]:])

    fw.close()    
    
def get_external_coord(filein):
    with open(filein, 'r') as f:
        [atoms, deriva, charge ,spin] = [int(s) for s in f.readline().split()]
        ele = np.zeros((atoms,), dtype = int)
        
        atomlist = np.zeros((atoms,3), dtype = float)

        for i in range(atoms):
            axestr = f.readline().split()
            ele[i] = int(axestr[0])
            atomlist[i][0] = float(axestr[1])*0.52917720859
            atomlist[i][1] = float(axestr[2])*0.52917720859
            atomlist[i][2] = float(axestr[3])*0.52917720859
    return ele, atomlist, deriva

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

if __name__ == '__main__':
    print(len(sys.argv))
    print(sys.argv)
    if len(sys.argv) == 8:
        gjff = sys.argv[1]
        filein = sys.argv[3]
        fileout = sys.argv[4]
        GaussianOniom(gjff,filein,fileout)
    else:
        print('Gaussian external interface call ONIOM in Gaussian ')
        print('')
        print('1 input: model.gjf; ')
        print('methods defined in model gjf')
        print('Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile')
        print('More infomation see Gaussian website')
        sys.exit()
    
    
    