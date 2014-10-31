## @package QCALC
#  QCALC PYTHON module to get infos from QCALC output file
#

#author:  L.Desorgher
#
#History: 
#---------
#      25/07/2014     Creation

import hepunit as unit
import numpy as np
from utilities import *
import os  

qcalc_dir="../qcalc_data/august2012/"

## Define the directory where the QCALC data are located 
#
#
def SetQCALCDirectory(dir):
    global qcalc_dir
    qcalc_dir=dir

## Get atomic or nuclear mass of a given  nucleus
#
#       
def GetMassFromQCALC1(Z,A,nuclear=False):
    qcalc_dic=ReadQCALCResults(Z,A)
    AMASS=qcalc_dic['AtomicMass(AMU)']*unit.amu_c2
    mass=AMASS
    if nuclear:
        electron_binding_energy=qcalc_dic['BindingEnergy/A']*A*unit.keV
        mass=AMASS+electron_binding_energy
    return mass

## Get atomic or nuclear mass of a given  nucleus
#
#  
def GetMassFromQCALC(Z,A,nuclear=False):
    #All the mass and excess mass in QCalc are for atomic mass
    #The binding energy gives the nuclear binding energy plus porbably the atomic electon binding energy 
    qcalc_dic=ReadQCALCResults(Z,A)
    AMASS=qcalc_dic['AtomicMass(AMU)']*unit.amu_c2
    mass=AMASS
    if nuclear:
        mass=Z*unit.proton_mass_c2+(A-Z)*unit.neutron_mass_c2-A*qcalc_dic['BindingEnergy/A']*unit.keV
        mass+=( 14.4381*np.power( Z , 2.39 ) + 1.55468*1e-6*np.power( Z , 5.35 ) )*unit.eV
        #print "MASS",mass/uni
    return mass


## Get Q value for alpha decay of a given  nucleus
#
#  
def GetAlphaQValueFromQCalc(Z,A):
    mass=GetMassFromQCALC(Z,A)
    massDaughter=GetMassFromQCALC(Z-2,A-4)
    print mass,massDaughter
    massAlpha=GetMassFromQCALC(2,4)
    if massDaughter>0.:
        return mass-massDaughter-massAlpha
    return None



## Get Q value for beta minus decay of a given  nucleus
#
# 
def GetBetaMinusQValueFromQCalc(Z,A):
    mass=GetMassFromQCALC(Z,A,nuclear=True)
    massDaughter=GetMassFromQCALC(Z+1,A,nuclear=True)
    print mass,massDaughter
    massElectron=unit.electron_mass_c2
    if massDaughter>0.:
        return mass-massDaughter-massElectron
    return None

## Get Q value for beta plus decay of a given  nucleus
#
# 
def GetBetaPlusQValueFromQCalc(Z,A):
    mass=GetMassFromQCALC(Z,A,nuclear=True)
    massDaughter=GetMassFromQCALC(Z-1,A,nuclear=True)
    print mass,massDaughter
    massElectron=unit.electron_mass_c2
    if massDaughter>0.:
        return mass-massDaughter-massElectron
    return None

## Read QCALC results
#
# 
def ReadQCALCResults(Z,A):
    file_name=qcalc_dir+"/"+"z%i.a%i" %(Z,A)
    file_name=file_name.replace("//","/")
    qcalc_dic={}
    if os.path.exists(file_name):
        qcalc_file_obj=open(file_name,'r')
        lines=qcalc_file_obj.readlines()
        for line in lines:
            words=line.split()
            if is_number(words[1]):
                qcalc_dic[words[0]]=float(words[1])
            if is_number(words[2]):
                qcalc_dic[words[0]+"_err"]=float(words[2])
    return qcalc_dic

if __name__ == "__main__":
    Z=92
    A=240
    print GetMassFromQCALC(Z,A,nuclear=False)
    print GetMassFromQCALC(Z,A,nuclear=True)

    