## @package Nuclear Wallet Cards
#  PYTHON module to read/ get info from nuclear wallet file
#

#author:  L.Desorgher
#
#History: 
#---------
#      21/08/2012     Creation

import hepunit as unit
import numpy as np
from utilities import *
import os
import QCalc

nuclear_wallet_file_name="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt"


## Set Nuclear Wallet File 
#  
# @param file_name path of the nuclear wallet cards 
#
def SetNuclearWalletFileName(file_name):
    global nuclear_wallet_file_name
    nuclear_wallet_file_name=file_name

## Get name of a given nucleus
#
#  @param Z atomic number
#  @param A mass number
def GetNucleusName(Z,A):
    map_AZ_to_names, map_names_to_AZ,A_vec,Z_vec = ReadNucleusNameList()
    if map_AZ_to_names is not None:
        if (A*1000+Z) in map_AZ_to_names:
            return map_AZ_to_names[A*1000+Z]
        else :
            return None

## Read the list of names of nucleus given in the nuclear wallet cards
#
#
def ReadNucleusNameList():
    file_object=open(nuclear_wallet_file_name,'r')
    lines=file_object.readlines()
    file_object.close()
    map_AZ_to_names={}
    map_names_to_AZ={}
    A_vec=[]
    Z_vec=[]
    
    for line in lines[1:]:
        words=line.split()
        A=int(words[0])
        Z=int(words[2])
        nr=A*1000+Z
        if nr not in map_AZ_to_names:
            name=words[0]+words[1].upper()
            map_AZ_to_names[nr]=name
            map_names_to_AZ[name]=[Z,A]
            A_vec+=[A]
            Z_vec+=[Z]
    return map_AZ_to_names,map_names_to_AZ,A_vec,Z_vec

## Read nuclear wallet file
#
#
def ReadNuclearWalletFile():
    file_object=open(nuclear_wallet_file_name,'r')
    lines=file_object.readlines()
    first_line=lines[0]
    multi_words_tags=["Mass Exc","T1/2 (txt)","T1/2 (seconds)","Dec Mode","Branching (%)"]
    new_first_line=first_line
    for tag in multi_words_tags:
        new_tag=tag.replace(" ","_")
        new_first_line=new_first_line.replace(tag,new_tag)
    indices=[]
    for tag in new_first_line.split():
        indices+=[new_first_line.find(tag)]
    indices+=[-1]
    tags=new_first_line.split()
    nuclear_info_dict={}
    for line in lines[1:]:
        for i in range(len(tags)):
            tag=tags[i]
            if tag not in nuclear_info_dict:
                nuclear_info_dict[tag]=[]
            nuclear_info_dict[tag]+=[line[indices[i]:indices[i+1]-1]]
    list_stable_nuclei=[]
    list_radioactive_nuclei=[]
    radioactive_nuc_dic={}
    stable_nuc_dic={}
    for i in range(len(nuclear_info_dict['A'])):
        A=int(nuclear_info_dict['A'][i])
        Z=int(nuclear_info_dict['Z'][i])
        Element=nuclear_info_dict['Element'][i]
        Energy_str=nuclear_info_dict['Energy'][i]
        Energy=string_to_float(nuclear_info_dict['Energy'][i],fill_value=-1.)
        Mass=A*unit.amu_c2+string_to_float(nuclear_info_dict["Mass_Exc"][i],fill_value=-999999999999999999999.)
        T_half=1e9
        T_half_str=nuclear_info_dict['T1/2_(seconds)'][i]
        is_stable=False
        name="z%i.a%i" %(Z,A)
        dic_to_fill=None
        if ("STABLE" in nuclear_info_dict['T1/2_(txt)'][i]):
            is_stable=True
            list_stable_nuclei+=[name]
            stable_nuc_dic[name]={'Z':Z,'A':A,'Element':Element,'Mass':Mass}
        else:
            T_half=string_to_float(T_half_str,fill_value=0.)
            Dec_Mode=nuclear_info_dict['Dec_Mode'][i].replace(" ","")
            b_ratio_str=nuclear_info_dict['Branching_(%)'][i]
            b_ratio=string_to_float(b_ratio_str.replace('<=',"").replace('<',"").replace('~',""),fill_value=-1.)
            if b_ratio_str.replace(" ","") =="":
                b_ratio=100.
            if name not in list_radioactive_nuclei:
                list_radioactive_nuclei+=[name]
                radioactive_nuc_dic[name]={'Z':Z,'A':A,'Element':Element,"excited_energy_vec":[],"Mass":Mass}
            decay=radioactive_nuc_dic[name]
            if Energy not in decay:
                decay[Energy]={}
                decay["excited_energy_vec"]+=[Energy]
            decay[Energy][Dec_Mode]=b_ratio
            decay[Energy]['T1/2']=T_half
    return list_stable_nuclei, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic

## Get atomic (default) or nuclear mass of an atom/nucleus
#
#   @param Z atomic number
#   @param A mass number
#   @param nuclear True/False  the nuclear/atomic mass is computed
def GetMassFromWalletCard(Z,A,nuclear=False):
    list_stable_nuclei, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic =ReadNuclearWalletFile()
    nuc_dic=None
    mass=-9999999999.
    
    if name in list_stable_nuclei:
         nuc_dic=stable_nuc_dic[name]
    if name in list_radioactive_nuclei:
         nuc_dic=radioactive_nuc_dic[name]
    if nuc_dic is not None:
        mass=nuc_dic['Mass']
        if nuclear:
            mass=mass-Z*unit.electron_mass_c2
            mass+=( 14.4381*np.power( Z , 2.39 ) + 1.55468*1e-6*np.power( Z , 5.35 ) )*unit.eV
    return mass

## Get Q value for alpha decays (atomic mass are considered) 
#
#   @param Z atomic number
#   @param A mass number
def GetAlphaQValueFromWalletCard(Z,A):
    mass=GetMassFromWalletCard(Z,A)
    massDaughter=GetMassFromWalletCard(Z-2,A-4)
    massAlpha=GetMassFromWalletCard(2,4)
    if massDaughter>0.:
        return mass-massDaughter-massAlpha
    return None

## Get Q value for beta minus decay (atomic mass are considered) 
#
#   @param Z atomic number
#   @param A mass number
def GetBetaMinusQValueFromWalletCard(Z,A):
    mass=GetMassFromWalletCard(Z,A)
    massDaughter=GetMassFromWalletCard(Z+1,A)
    massElectron=unit.electron_mass_c2
    if massDaughter>0.:
        return mass-massDaughter-massElectron
    return None

## Get Q value for beta minus decay (atomic mass are considered) 
#
#   @param Z atomic number
#   @param A mass number
def GetBetaPlusQValueFromWalletCard(Z,A):
    mass=GetMassFromWalletCard(Z,A)
    massDaughter=GetMassFromWalletCard(Z-1,A)
    massElectron=unit.electron_mass_c2
    if massDaughter>0.:
        return mass-massDaughter-massElectron
    return None

## Read radioactive decay data  of a given nucleus
#
#   @param Z atomic number
#   @param A mass number     
def ReadDecaysInWalletCards(Z,A):    
    
    list_stable_nuclei, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic =ReadNuclearWalletFile()
    name = "z%i.a%i" %(Z,A)
    Parent_elevel_list=[]
    Parent_half_life_list =[]
    allowed_decays=[]
    decay_vec=[]
    qcalc_dic=QCalc.ReadQCALCResults(Z,A)
    q_map={}
    if 'Qalpha' in qcalc_dic:
        q_map['A']=qcalc_dic['Qalpha']

            
    if 'Qbeta-' in qcalc_dic:
        q_map['B-']=qcalc_dic['Qbeta-']
    
    if 'QEC' in qcalc_dic:
        q_map['EC']=qcalc_dic['QEC']
    
    if 'Qbeta+' in qcalc_dic:
        q_map['Qbeta+']=qcalc_dic['Qbeta+']
        #q_map['EC']=q_map['Qbeta+']
    q_map['IT']=0.
    if name in list_radioactive_nuclei:
        nucleus_dic=radioactive_nuc_dic[name]
        map={"A":"Alpha","B-":"BetaMinus","EC":"BetaPlus","IT":"IT"}
        Parent_elevel_list=nucleus_dic['excited_energy_vec'] 
        for elevel in Parent_elevel_list:
            decays=nucleus_dic[elevel]
            Parent_half_life_list+=[decays["T1/2"]]
            decay={'type_list':[],'Q_list':[],'BR_list':[]}
            decay_valid=False
            for dec_type in map :
                if dec_type in decays and dec_type in q_map:
                    br=decays[dec_type]
                    dec_name=map[dec_type]
                    q=q_map[dec_type]+elevel    
                       
                    #We only consider observed decay
                    ################################
                    if br>0. and q>0.:
                        decay_valid=True
                        if dec_name not in allowed_decays:
                            allowed_decays+=[dec_name]
                        decay['type_list']+=[dec_name]
                        decay['BR_list']+=[decays[dec_type]]
                        q=q_map[dec_type]+elevel    
                        decay['Q_list']+=[q]
            decay_vec+=[decay]
                    
    return allowed_decays,decay_vec,Parent_elevel_list,Parent_half_life_list
 

if __name__ == "__main__":
    Z=91
    A=239
    print ReadDecaysInWalletCards(Z,A)
         