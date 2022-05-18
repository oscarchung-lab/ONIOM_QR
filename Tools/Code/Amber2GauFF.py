#!/usr/bin/env python3
import sys
import os
import numpy as np

def amberff2prm(ffname, prmname):
    #split ff file to atom, bond, angle, torsion, dihe and vdw part
    os.system("grep '^[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]*' "+ffname+" > atom.tmp")
    os.system("grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' "+ffname+" > bond.tmp")
    os.system("grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' "+ffname+" > angle.tmp")
    os.system("grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' "+ffname+" > torsion.tmp")
    os.system("grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' "+ffname+" > improper.tmp")
    os.system("grep '^  [A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' "+ffname+" > vdw.tmp")
    #preprare write the prm file
    fw=open(prmname, 'w')
    fw.write('! \n')
    fw.write('! Non-bonded interaction \n')
    fw.write('! \n')
    fw.write('NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000 \n')
    fw.write('! \n')
    fw.write('! Stretches \n')
    fw.write('! \n')
    numbond=0
    with open('bond.tmp', 'r') as f:
        while True:
            bondline=f.readline()
            if not bondline:
                break 
            bond1=bondline[0:12].replace("-"," ").rstrip()
            bond2=bondline[13:23].rstrip().lstrip()
            fw.write('HrmStr1 '+bond1+'  '+bond2+'\n')
            numbond=numbond+1
    print('number of bond: ', numbond)
    fw.write('! \n')
    fw.write('! Angles \n')
    fw.write('! \n')

    numangle=0
    with open('angle.tmp', 'r') as f:
        while True:
            angleline=f.readline()
            if not angleline:
                break
            angle1=angleline[0:17].replace("-"," ").rstrip().lstrip()
            angle2=float(angleline[18:30].rstrip().lstrip())
            fw.write('%24s%9.2f\n' % ('HrmBnd1 '+angle1, angle2))
            numangle=numangle+1
    print('number of angle: ', numangle)
    fw.write('! \n')
    fw.write('! Torsions \n')
    fw.write('! \n')

    numtorsion=0
    with open('torsion.tmp', 'r') as f:
        tortmp=' '
        torpara1=np.zeros(4, dtype = np.float32)
        torphase1=np.zeros(4, dtype = np.int32)
        while True:
            torline=f.readline()
            if not torline:
                if tortmp!=' ':
                    toreles = tortmp.replace("X ","* ").replace("-"," ")
                    fw.write('%18s%4d%4d%4d%4d%7.3f%7.3f%7.3f%7.3f%4.1f\n' %('AmbTrs '+toreles, torphase1[0], torphase1[1], torphase1[2], torphase1[3], torpara1[0], torpara1[1], torpara1[2], torpara1[3], divider))
                    numtorsion=numtorsion+1
                break
            if torline[0:1] == "X":
                if tortmp!=' ':
                    toreles = tortmp.replace("X ","* ").replace("-"," ")
                    fw.write('%18s%4d%4d%4d%4d%7.3f%7.3f%7.3f%7.3f%4.1f\n' %('AmbTrs '+toreles, torphase1[0], torphase1[1], torphase1[2], torphase1[3], torpara1[0], torpara1[1], torpara1[2], torpara1[3], divider))
                    tortmp=' '
                    numtorsion=numtorsion+1
                toreles = torline[0:11].replace("X ","* ").replace("-"," ")
                axestr = torline[13:55].split()
                divider= float(axestr[0])
                #barterm= float(axestr[1])
                #phase= float(axestr[2])
                perio= int(float(axestr[3]))
                torpara=np.zeros(4, dtype = np.float32)
                torphase=np.zeros(4, dtype = np.int32)
                torpara[perio-1]=float(axestr[1])
                torphase[perio-1]=int(float(axestr[2])+0.5)
                fw.write('%18s%4d%4d%4d%4d%7.3f%7.3f%7.3f%7.3f%4.1f\n' %('AmbTrs '+toreles, torphase[0], torphase[1], torphase[2], torphase[3], torpara[0], torpara[1], torpara[2], torpara[3], divider))
                numtorsion=numtorsion+1
            else:
                if tortmp==' ':
                    tortmp = torline[0:11]
                    axestr = torline[12:55].split()
                    divider= float(axestr[0])
                    
                    perio= int(abs(float(axestr[3])))
                    torpara1[perio-1]=float(axestr[1])
                    torphase1[perio-1]=int(float(axestr[2])+0.5)
                elif torline[0:11] == tortmp:
                    axestr = torline[12:55].split()
                    #divider= float(axestr[0])

                    perio= int(abs(float(axestr[3])))
                    torpara1[perio-1]=float(axestr[1])
                    torphase1[perio-1]=int(float(axestr[2])+0.5)
                else:
                    toreles = tortmp.replace("X ","* ").replace("-"," ")
                    fw.write('%18s%4d%4d%4d%4d%7.3f%7.3f%7.3f%7.3f%4.1f\n' %('AmbTrs '+toreles, torphase1[0], torphase1[1], torphase1[2], torphase1[3], torpara1[0], torpara1[1], torpara1[2], torpara1[3], divider))
                    numtorsion=numtorsion+1
                    tortmp = torline[0:11]
                    torpara1[:]=0.0
                    torphase1[:]=0
                    axestr = torline[12:55].split()
                    divider= float(axestr[0])
                    
                    perio= int(abs(float(axestr[3])))
                    torpara1[perio-1]=float(axestr[1])
                    torphase1[perio-1]=int(float(axestr[2])+0.5)
    print('number of torsion: ', numtorsion)          
    fw.write('! \n')
    fw.write('! Improper torsions \n')
    fw.write('! \n')            

    numimproper=0
    with open('improper.tmp', 'r') as f:
        while True:
            improperline=f.readline()
            if not improperline:
                break
            improper1=improperline[0:11].replace("X ","* ").replace("-"," ")
            axestr = improperline[13:55].split()
            improper2= float(axestr[0])
            improper3= float(axestr[1])
            improper4= float(axestr[2])
            fw.write('%18s%6.1f%7.1f%4.1f\n' % ('ImpTrs '+improper1, improper2,improper3,improper4))
            numimproper=numimproper+1
    print('number of improper: ', numimproper)  
    fw.write('! \n')
    fw.write('! Vanderwaals parameters \n')
    fw.write('! \n')            

    numvdw=0
    with open('vdw.tmp', 'r') as f:
        while True:
            vdwline=f.readline()
            if not vdwline:
                break
            vdw1=vdwline[2:4]
            vdw2=vdwline[14:20]
            vdw3=vdwline[22:34]
            
            fw.write('VDW '+vdw1+' '+vdw2+'  '+vdw3+'\n')
            numvdw=numvdw+1
    print('number of vdw: ', numvdw) 
    os.system("rm -f *.tmp")

if __name__ == '__main__':
    print(len(sys.argv))
    if len(sys.argv) == 2:
        
        amberff = sys.argv[1]
        fname=os.path.splitext(amberff)[0]
        gauprm=fname+'.prm'
        amberff2prm(amberff,gauprm)
    else:
        print('Inputs: Amber FF parameter file (.dat)')


