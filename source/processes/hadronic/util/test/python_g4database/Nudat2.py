#
#ENSDF python module to read,interprets ENSDF, and produce G4 raddecay and evaporation database  file
#
#author:  L.Desorgher
#
#History: 
#---------
#      21/08/2012     Creation

import hepunit as unit
import numpy 
import zipfile
import os
import subprocess
import numpy as np
import numpy

from QCalc import *
from NuclearWalletCards import *
from utilities import *


uma_in_MeV=931.494028
nudat2_dir="../nudat2_data/nuclear_decay/august2012"

## Define the directory where nudat2 files are located
#
#  @param dir directory where the nudat2 data are located
# 
def SetDataDirectory(dir):
    global nudat2_dir
    nudat2_dir=dir
    

def GetRadiationIntensitiesFromRadDecayNudat2(Z,A):
    global ic,Avec,Zvec,HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec
    ic=0
    file_list=np.array(os.listdir(nudat2_dir))
    Z_vec=[]
    for file_name in file_list:
        Z_vec+=[int(file_name.split("[")[1].split(",")[0])]
    Z_vec=np.array(Z_vec)
    indices=np.argsort(Z_vec)
    Z_vec=Z_vec[indices]
    file_vec=file_list[indices]
    ind=np.where(Z>=Z_vec)[0][-1]
    Avec,Zvec,HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec=np.loadtxt("%s/%s" %(dir,file_vec[ind]), delimiter='\t',skiprows=1, usecols=(0,2,6,8,9,10,14),unpack=True,
                                                                           converters={8:type_radiation_to_float,9:sub_type_radiation_to_float,10:string_to_float,14:string_to_float})
    ind=np.where((Avec==A)*(Zvec==Z))
    return HalfTime_vec[ind],TypeRad_vec[ind],SubTypeRad_vec[ind],Energy_vec[ind],Intensity_vec[ind]
    

def GetGammaIntensitiesFromRadDecayNudat2(Z,A,with_xrays=False):
    import Nudat2WebGetter as nudat2
    theNudat2Getter= nudat2.Nudat2WebDataGetter()
    HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec=GetRadiationIntensitiesFromRadDecayNudat2(Z,A)
    print HalfTime_vec
    if len(HalfTime_vec)<=0:
        return np.array([]),np.array([])
    tmax=max(HalfTime_vec)
        
    ind=np.where((TypeRad_vec == 5) * (SubTypeRad_vec ==-1.) * (HalfTime_vec == tmax))
    return Energy_vec[ind],Intensity_vec[ind]

def GetAlphaIntensitiesFromRadDecayNudat2(Z,A):
    HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec=GetRadiationIntensitiesFromRadDecayNudat2(Z,A)
    ind=np.where((TypeRad_vec == 0) )
    return Energy_vec[ind],Intensity_vec[ind]
    
def GetXrayIntensitiesFromRadDecayNudat2(Z,A):
    HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec=self.GetRadiationIntensitiesFromRadDecayNudat2(Z,A)
    ind=np.where((TypeRad_vec == 5) * (SubTypeRad_vec ==0.))
    return Energy_vec[ind],Intensity_vec[ind]
    
def GetConversionElectronIntensitiesFromRadDecayNudat2(Z,A):
    HalfTime_vec,TypeRad_vec,SubTypeRad_vec,Energy_vec,Intensity_vec=self.GetRadiationIntensitiesFromRadDecayNudat2(Z,A)
    ind=np.where((TypeRad_vec == 6) * (SubTypeRad_vec ==2))
    return Energy_vec[ind],Intensity_vec[ind]
    
            
            
    
    
            
        
        
        
        

    
    
    






                 
   
            
            
        

    
     