## @package utilities
#  Some utilities for the production of g4database
#

#author:  L.Desorgher
#
#History: 
#---------
#      25/07/2014     Creation

import hepunit as unit
import numpy as np


## Convert nudat2 type of radiation in a float
#
#
def type_radiation_to_float(type_rad):
    #print type_rad
    string_rad=type_rad.replace(" ","")
    dic= {'A':0.,'BM':1.,'BMAV':2.,'BP':3.,'BPAV':4.,'G':5.,'E':6.}
    if string_rad in dic:
        return dic[string_rad]
    else:
        return -1.


## Convert nudat2 sub type of radiation in a float
#
#
def sub_type_radiation_to_float(sub_type_rad):
    global ic
    ic+=1
    string_sub=sub_type_rad
    words=sub_type_rad.split()
    if len(words)>0:
        string_sub=words[0]
    string_sub=string_sub.replace(" ","")
    dic= {'XR':0.,'Auger':1.,'CE':2.}
    if string_sub in dic:
        return dic[string_sub]
    else:
        return -1.

## Convert K,L1,..,M1.. shell label in an integer 
#
#
def shell_label_to_id(label):
    shell_label=label.replace(" ","")
    shells_id={"K":0,"L1":1,"L2":2,"L3":3,"M1":4,"M2":5,"M3":6,"M4":7,'M5':8}
    if (shell_label in shells_id):
        return shells_id[shell_label]
    elif shell_label.find("-tot") !=-1:
        return -1
    elif shell_label == "IPF":
        return 10
    elif shell_label == "Tot":
        return -1
    elif shell_label.find("/") !=-1 :
        return -1
    else : #Other shells
        return 9
## Check if  the input is a number
#
#
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

## Convert a string to a float
#
#
def string_to_float(s,fill_value=0.):
    s1=s.replace("**********","0").replace(">","")
    try:
        r=float(s1)
        return r
    except ValueError:
        return fill_value

## View the contain of a file 
#
#   
def ViewFile(file_name):
    print file_name
    if not os.path.exists(file_name):
        print "File does not exist"
        return
    file_object=open(file_name,'r')
    lines=file_object.readlines()
    for line in lines:
        print line[:-1]


## Compute the change in Z and A for a given decay type
#
# 
def GetDeltaZandAForDecayType(decay_type):
    if  decay_type in  ['MshellEC',"KshellEC","BetaPlus","LshellEC"]:
        return -1,0
    if  decay_type in  ['BetaMinus']:
        return 1,0
    if  decay_type in  ['Alpha']:
        return -2,-4
    if  decay_type in  ['IT']:
        return 0,0

## Compute the kinetic energy of the daughters of a two body decay
#
#
def ComputeEnergyProductsTwoBodyDecay(parentMass,daughter1Mass,daughter2Mass):
    Q=parentMass-daughter1Mass-daughter2Mass
    daughter1Ekin=(daughter2Mass*Q+Q*Q/2.)/(parentMass)
    daughter2Ekin=(daughter1Mass*Q+Q*Q/2.)/(parentMass)
    return daughter1Ekin,daughter2Ekin     


## Compute delta JPi and possible multi-polarity for a gamma decay
#
#
def GetDJpiAndPossibleMultipolarity(JPi_vec1,JPi_vec2):
    DJpi_str=None
    MultiVec=None
    J1=None
    Pi1=None
    J2=None
    Pi2=None
    if (JPi_vec1 is not None and len(JPi_vec1)>0):
        J1=JPi_vec1[0][0]
        if len(JPi_vec1[0])>1:
            Pi1=JPi_vec1[0][1]
    if (JPi_vec2 is not None and len(JPi_vec2)>0):
        J2=JPi_vec2[0][0]
        if len(JPi_vec2[0])>1:
            Pi2=JPi_vec2[0][1]
        
        
        
    if (J1 is not None and J2 is not None and Pi1 is not None and Pi2 is not None):
        DJpi_str="%i+" %(abs(J1-J2))
        if (Pi1*Pi2) <0:
            DJpi_str=DJpi_str.replace("+","-")
        if (J1>0 and DJpi_str=="0-"):
            DJpi_str="1-"
        if (J1>0 and DJpi_str=="0+"):
            DJpi_str="1+"
        
        
        lmin=abs(J1-J2)
        lmax=J1+J2
        if lmax>=1 and lmin<=0:
             lmin=1
        L_vec=[]
        if (lmax >=lmin):
            L_vec=lmin+np.arange(lmax-lmin+1)
        
        MultiVec=[]
        for L in L_vec:
            MultiType="E"
            if (L%2): #odd case
                MultiType="M"
                if (Pi1*Pi2) <0:
                    MultiType="E"
            else:
                MultiType="E"
                if (Pi1*Pi2) <0:
                    MultiType="M"
                
            MultiVec+=["%s%i" %(MultiType,L)]
    return  DJpi_str, MultiVec
