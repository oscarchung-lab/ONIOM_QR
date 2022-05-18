#!/usr/bin/env python3
import sys
import os
import numpy as np
from string import digits

def mol2gjf(mol2fname, gjffname):
    fw = open(gjffname,"w")

    with open(mol2fname, 'r') as f:
        f.readline()
        f.readline()
        axestr = f.readline().split()
        numa=int(axestr[0])

        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        table = str.maketrans('', '', digits)
        for i in range(numa):
            axestr = f.readline().split()
            ele = axestr[1].translate(table)
            x = float(axestr[2])
            y = float(axestr[3])
            z = float(axestr[4])
            atomt = axestr[5].upper()
            atomchg = axestr[8]

            fw.write('%-18s%1s%14.8f%14.8f%14.8f%2s\n' % (' '+ele+'-'+atomt+'-'+atomchg,'0',x,y,z,' H'))


if __name__ == '__main__':
    print(len(sys.argv))
    if len(sys.argv) == 2:
        
        mol2f = sys.argv[1]
        fname=os.path.splitext(mol2f)[0]
       
        mol2gjf(mol2f, fname+'.dat')
    else:
        print('Inputs: RESP mol2 file from anterchamber (.mol2)')


