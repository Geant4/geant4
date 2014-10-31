## @package ENSDF
#  ENSDF PYTHON module to read and interpret ENSDF files
#
#  It reads ENSDF files and get the necessary information for building the G4Radiaoctive and G4PhotoEvaporation database


#author:  L.Desorgher
#
#History: 
#---------
#      21/08/2014     Creation

import hepunit as unit
import numpy 
import zipfile
import os
import subprocess
import numpy as np
import numpy

import QCalc
import NuclearWalletCards 
from utilities import *


global adopted_levels_and_gammas,ids_isomer_to_add_in_database,elevels_isomer_to_add_in_database,elevels_daughter_isomer_to_add_in_database,half_life_isomer_to_add_in_database 
adopted_levels_and_gammas={}
ids_isomer_to_add_in_database=[]
elevels_isomer_to_add_in_database=[]
elevels_daughter_isomer_to_1add_in_database=[]
half_life_isomer_to_add_in_database=[]
uma_in_MeV=931.494028
zip_file_name_decays="../ensdf_data/september_2013/decays.zip"
zip_file_name_levels="../ensdf_data/november_2013/all_levels.zip"
zip_file_name_decays_xundl="../ensdf_data/xundl.zip"



## Set the decay and level zip files that contain the ENSDF data
#
#  @param path_decay_file zip file containing the whole ENSDF files for decays
#  @param zip_file_name_levels zip file containing the whole ENSDF files for levels  
def SetZipFilesENSDF(path_decay_file,path_level_file):
    global zip_file_name_decays,zip_file_name_levels
    zip_file_name_decays=path_decay_file
    zip_file_name_levels=path_level_file

## Set the decay zip file that contain the XUNDL data
#
#  @param path_decay_file zip file containing the whole ENSDF files for decays 
def SetZipFileXUNDL(path_decay_file):
    global zip_file_name_decays_xundl,zip_file_name_levels_xundl
    zip_file_name_decays_xundl=path_decay_file


## Read the list of available nuclei in ENSDF 
#
def ReadListOfAvailableNucleiInENDSF():
        theZipFile=zipfile.ZipFile(zip_file_name_decays)
        list_decay_names=theZipFile.namelist()
        theZipFile.close()
        
        theZipFile=zipfile.ZipFile(zip_file_name_levels)
        list_level_names=theZipFile.namelist()
        theZipFile.close()
        return list_decay_names,list_level_names


#Function to read values of fields
#############################

## Get value of Energy field in ENSDF record 
#
#
def GetValueForEField(efield_str,de_field_str=None):
    aString=efield_str
    directly_measured=True
    if "(" in aString:
        directly_measured=False
        aString=aString.replace("(","").replace(")","")
        
    is_given=False
    value=0.
    value_type="Unknown"
    ref_level="Zero"
    if  is_number(aString):
        is_given=True
        value=float(aString)
        value_type="Absolute"
    elif "+" in aString:
        words=aString.split("+")
        if words[0] in ["SP","SN"] and is_number(words[1]):
            value=float(words[1])
            is_given=True
            value_type="Resonance"
        elif is_number(words[0]):
            is_given=True
            value=float(words[0])
            value_type="AboveRefValue_"+words[1]
            ref_level=words[1].replace(" ","")
        elif is_number(words[1]):
            is_given=True
            value=float(words[1])
            value_type="AboveRefValue_"+words[0]
            ref_level=words[0].replace(" ","")
    else:
        ref_level=aString.replace(" ","")
        if aString in ["X","Y","Z"]:
            value_type="AboveRefValue_"+aString
            value=0.
            
    #value*=unit.keV
    value_unit="keV"
    if de_field_str is None:
        return value,is_given,value_type,value_unit,directly_measured,ref_level
    else:
        de_value,de_value_given,is_upper_limit,is_lower_limit=GetValueForUncertaintyField(de_field_str)
        return value,is_given,value_type,value_unit,directly_measured,ref_level, \
                            de_value,de_value_given,is_upper_limit,is_lower_limit

## Get value of uncertainty field in ENSDF record 
#
#       
def GetValueForUncertaintyField(field):
    is_upper_limit=False
    is_lower_limit=False
    value_given=False
    value=0.
    if is_number(field):
        value=float(field)
        value_given=True
    elif field in ["LE","LT"]:
        is_upper_limit=True
    elif field in ["GE","GT"]:
        is_lower_limit=True
    return  value,value_given,is_upper_limit,is_lower_limit

## Get JPi from string 
#
#          
def GetJPiFromString(Jpi_str):#Jpi string are either J, pi, orJpi
    J=None
    Pi=None
    #print "JPi",Jpi_str
    number_str=Jpi_str.replace("+"," ").replace("-"," ").replace("/","  ")
    words=number_str.split()
    #print words
    if len(words)>1:
        number_str=words[0]
        
    if is_number(number_str):
        J=float(number_str)
        if ("/" in Jpi_str):
            J=J/2.
            
    if "+" in Jpi_str:
        Pi=1
    if "-" in Jpi_str:
        Pi=-1
    return J,Pi

## Get possible Jpi values from Jpi field in ENSDF record 
#
# 
def GetValuesForJPiField(field): #JPi
    JPi_vec=[]
    ref_level=None
    field_str=field.replace("(","").replace(")","").replace("]","").replace("[","")

    if "OR" in field_str or "," in field_str or "AND" in field_str:
        vec_str=field_str.replace("OR"," ").replace(","," ").replace("AND","").split()
        #print vec_str
        JPi_vec=[]
        for str1 in vec_str:
            J,Pi=GetJPiFromString(str1)
            if J is not None:
                JPi_vec+=[[J]]
            if Pi is not None and len(JPi_vec)>0:
                JPi_vec[-1]+=[Pi]
    elif "NOT" in field_str or "NATURAL" in field_str:
        return [],None
    elif "AP" in field_str:
        J,Pi=GetJPiFromString(field_str.replace("AP"," "))
        if J is not None:
            JPi_vec=[[J]]
            if Pi is not None:
                JPi_vec=[[J,Pi]]
    elif "LE" in field_str:
        J,Pi=GetJPiFromString(field_str.replace("AP"," "))
        if J is not None:
            while(J>=0):
                JPi_vec+=[[J]]
                if Pi is not None:
                    JPi_vec[-1]=[J,Pi]
                J=J-1
    elif "GE" in field_str:
        J,Pi=GetJPiFromString(field_str.replace("AP"," "))
        if J is not None:
            for i in range(3):
                JPi_vec+=[[i+J]]
                if Pi is not None:
                    JPi_vec[-1]=[J,Pi]
    elif "TO" in field_str or ":" in field_str:
        JPi_strings=field_str.replace("TO"," ").replace(":"," ").split()
        #print JPi_strings
        J1,Pi1=GetJPiFromString(JPi_strings[0])
        #print J1,Pi1
        J2,Pi2=GetJPiFromString(JPi_strings[1])
        #print J2,Pi2
        J=J1
        JPi_vec=[]
        while(J<=J2):
            JPi_vec+=[[J]]
            J+=1
        if (Pi1 is None ):
            if (Pi2 is not None):
                for i in range(len(JPi_vec)):
                    JPi_vec[i]+=[Pi2]
        else :
            JPi_vec[0]+=[Pi1]
            for  i in range(len(JPi_vec)-1):
                    JPi_vec[i+1]+=[1,-1]
            if (Pi2 is not None):
                JPi_vec[-1]=[JPi_vec[-1][0],Pi2]
    else:
        level_found=False
        i=0
        ref_vec=["K","L","M","N","O","P","Q"]
        JPi_str=field_str
        while not level_found and i<len(ref_vec):
            test_str=ref_vec[i]+"+"
            if test_str in field_str:
                level_found=True
                ref_level=ref_vec[i]
                JPi_str=JPi_str.replace(test_str,"")
            i+=1
        J,Pi=GetJPiFromString(JPi_str)
        if J is not None:
            JPi_vec=[[J]]
            if Pi is not None:
                JPi_vec=[[J,Pi]]
    return JPi_vec,ref_level

## Get half time value  from T field in ENSDF record 
#
#      
def GetValueForTField(tfield_str): 
    half_time=0.
    is_stable=False
    determined=False
    words =tfield_str.split()
    if ("STABLE" in tfield_str):
        value=-99999999999
        determined=True
        is_stable=True
    if len(words)>1:
        determined=True
        half_time=float(words[0])
        unit_str=words[1]
        unit_value=0.
        hour=3600.
        minute=60.
        day=hour*24.
        year=365.25*day
        time_unit_dic={'MS':1e-3,'US':1e-6,'NS':1.e-9,'FS':1.e-15,'AS':1e-18,'PS':1e-12,'S':1,"Y":year,"D":day,"H":hour,"M":minute}
        energy_unit_dic={"MEV":unit.MeV,"KEV":unit.keV,"EV":unit.eV,'GEV':unit.GeV}
        if unit_str in time_unit_dic.keys():
            unit_value=time_unit_dic[unit_str]
        elif unit_str in energy_unit_dic:
            unit_value=4.55706*1.e-25*unit.GeV/energy_unit_dic[unit_str]
            half_time=1./half_time
        
        half_time*=unit_value
        
    return half_time,is_stable,determined

## Get beta decay type value from UN field in ENSDF record 
#
#    
def GetValueForUNField(unfield_str):#allowed , forbidden,....
    is_forbidden_decay=False
    is_forbidden_unique_decay=False
    forbidness_order=0
    if unfield_str !="  ":
        is_forbidden_decay=True
        forbidness_order=int(unfield_str[0])
        if (unfield_str[1]=="U"):
            is_forbidden_unique_decay=True
    g4_beta_type_flag=""
    if not is_forbidden_decay:
        g4_beta_type_flag="allowed"
    elif not is_forbidden_unique_decay and forbidness_order<=3:
        g4_beta_type_flag=["first","second","third"][forbidness_order-1]+"Forbidden"
    elif is_forbidden_unique_decay and forbidness_order<=3:
        g4_beta_type_flag="unique"+["First","Second","Third"][forbidness_order-1]+"Forbidden"
    return is_forbidden_decay,is_forbidden_unique_decay,forbidness_order,g4_beta_type_flag


## Get value from Number field in ENSDF record 
#
#  Different type of values can be contained in this field: QValue, Logft, ...
def GetValueForNumberField(number_field_str):#QValue,Logft,....
    value=0.
    value_is_given=False
    if (is_number(number_field_str)):
        value=float(number_field_str)
        value_is_given=True
    return value,value_is_given

## Get value from intensity field in ENSDF record 
#
#  Different types of intensity can be contained in this field: IA, IB, IE 
def GetValueForIntensityField(intensity_field_str):#IA,IB,IE....
    aString=intensity_field_str
    directly_measured=True
    if "(" in aString:
        directly_measured=False
        aString=aString.replace("(","").replace(")","")
    value=0.
    value_is_given=False
    if (is_number(aString)):
        value=float(aString)
        value_is_given=True
    return value,value_is_given,directly_measured

#         
#Methods for Reading of ENSDF records
#####################################


## Read parent record
#
#  
def ReadParentRecord(line):#no comment line with P in the 7th line
    nuc_name=line[0:5]
    e_level,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    Jpi_str=line[21:30].replace(" ","")
    half_life,is_stable,aBool=GetValueForTField(line[39:49])
    half_life_str=line[39:49]
    Qval,is_Qval_given=GetValueForNumberField(line[64:74])
    return nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given 

## Read normalization record
#
# 
def ReadNormalizationRecord(line):#no comment line with P in the 7th l
    nr,is_nr_given=GetValueForNumberField(line[9:19])
    nt,is_nt_given=GetValueForNumberField(line[21:29])
    br,is_br_given=GetValueForNumberField(line[31:39])
    nb,is_nb_given=GetValueForNumberField(line[42:49])
    np,is_np_given=GetValueForNumberField(line[55:62])
    return nr,nt,br,nb,np,is_nr_given,is_nt_given,is_br_given,is_nb_given,is_np_given

## Read alpha record
#
#   
def ReadAlphaRecord(line):
    energy,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    intensity,intensity_given,aBool=GetValueForIntensityField(line[21:29])
    
    hf,hf_given=GetValueForNumberField(line[31:29])
        
    return energy,e_value_unit,intensity,hf,e_is_given,intensity_given,hf_given

## Read level record
#
#
def ReadLevelRecord(line):
    energy,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    """
    if  "AboveRefValue_" in e_value_type:
        e_is_given=False
    """
    JPi_vec,x=GetValuesForJPiField(line[21:39])
    half_life,is_stable,aBool=GetValueForTField(line[39:49])
    meta_str=line[77:79]          
    return energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,line[21:39]
 
## Read beta minus record
#
#  
def ReadBetaMinusRecord(line):#no comment line with P in the 7th line
    end_energy,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    intensity,is_intensity_given,aBool=GetValueForIntensityField(line[21:29])
    logft,is_logft_given=GetValueForNumberField(line[41:49])
    is_forbidden_decay,is_forbidden_unique_decay,forbidness_order,g4_beta_type_flag=GetValueForUNField(line[77:79])
    return end_energy, e_value_unit,e_is_given,intensity,is_intensity_given,\
                    logft,is_logft_given,is_forbidden_decay,is_forbidden_unique_decay,forbidness_order,\
                    g4_beta_type_flag

## Read electron conversion (EC) record
#
#  
def ReadECRecord(line):#EC and Beta+ info
    energy,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    IB,is_IB_given,aBool=GetValueForIntensityField(line[21:29])   
    IE,is_IE_given,aBool=GetValueForIntensityField(line[31:39])  
    
    TI,is_TI_given,aBool=GetValueForIntensityField(line[64:74])  
    logft,is_logft_given=GetValueForNumberField(line[41:49])
    is_forbidden_decay,is_forbidden_unique_decay,forbidness_order,g4_beta_type_flag=GetValueForUNField(line[77:79])
              
    return energy, e_value_unit,e_is_given,IB,is_IB_given, \
            IE,is_IE_given,TI,is_TI_given,logft,is_logft_given, \
            is_forbidden_decay,is_forbidden_unique_decay,forbidness_order,g4_beta_type_flag


## Read gamma record
#
#  A gamma record defines a gamma emission from an excited state to a lower excited state observed in 
#  a given radioactive decay 
def ReadGammaRecord(line):
    energy,e_is_given,e_value_type,e_value_unit,e_directly_measured,ref_level=GetValueForEField(line[9:19])
    
    RI,is_RI_given,aBool=GetValueForIntensityField(line[21:29])
    MP_str=line[31:41].replace(" ","")
    
    MR,is_MR_given=GetValueForNumberField(line[41:49]) #mixing ratio 
    CC,is_CC_given=GetValueForNumberField(line[55:62]) #conversion coefficient
    
    TI,is_TI_given,aBool=GetValueForIntensityField(line[64:74])  
              
    return energy, e_value_unit,e_is_given,ref_level,RI,is_RI_given, \
            MP_str,MR,is_MR_given,CC, is_CC_given,TI,is_TI_given


## Read identification record
#
#  The identification record defines for a given ENSDF block the type of decay, if this 
#  block defines adopted levels of an isotope, or 
#  if the given information in the block is a tentative  
def ReadIdentificationRecord(line):
    nuc_id=line[0:5]
    dsid=line[9:39]
    is_beta_minus_decay=False
    is_EC_beta_plus_decay=False
    is_IT_decay=False
    is_alpha_decay=False
    is_tentative=False
    is_adopted_levels=False
    is_not_observed=False
    is_gammas=False
    if (" B- DECAY" in dsid):
        is_beta_minus_decay=True
    
    if (" A DECAY" in dsid):
        is_alpha_decay=True
    
    if (" EC DECAY" in dsid or " B+ DECAY" in dsid):
        is_EC_beta_plus_decay=True
    
    if (" IT DECAY" in dsid):
        is_IT_decay=True
    
    if ("TENTATIVE" in dsid):
        is_tentative=True
    if ("ADOPTED LEVELS" in dsid):
        is_adopted_levels=True
        if ("GAMMAS" in dsid):
            is_gammas=True
        if ("NOT OBSERVED " in dsid):
            is_not_observed=True
    
    
    return nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
          is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed  

## Read continuation record
#
#      
def ReadContinuationCard(line):

    words=line[9:].split("$")
    op_lists=["<=",">=","=","<",">"," EQ "," AP "," LT "," LE "," GT "," GE "]
    card_dic={}
    for word in words:
        used_ops=[]
        pos_ops=[]
        for op in op_lists:
            op_pos=word.find(op)
            if op_pos>=0:
                used_ops+=[op]
                pos_ops+=[op_pos]
        if len(used_ops)>0:
            used_ops=numpy.array(used_ops)
            pos_ops=numpy.array(pos_ops)
            indices=numpy.argsort(pos_ops)
            used_ops=used_ops[indices]
            pos_ops=pos_ops[indices]
            
            var_tag=word[0:pos_ops[0]].replace(" ","")
            card_dic[var_tag]={}
            for i in range(len(used_ops)):
                start_pos=pos_ops[i]+len(used_ops[i])
                value_str=word[start_pos:]
                if (i<len(used_ops)-1):
                    end_pos=pos_ops[i+1]+len(used_ops[i+1])
                    value_str=word[start_pos:end_pos]
                if (len(value_str.split())>0):
                    card_dic[var_tag]["value_str"]=value_str.split()[0]
                if (len(value_str.split())>0 and is_number(value_str.split()[0])):
                    val=float(value_str.split()[0])
                    if (used_ops[i] in ["="," EQ "," AP "]):
                        card_dic[var_tag]["value"]=val
                    elif (used_ops[i] in ["<=","<"," LE "," LT "]):
                        card_dic[var_tag]["max"]=val
                    elif (used_ops[i] in [">=",">"," GE "," GT "]):
                        card_dic[var_tag]["min"]=val
        
            
    return card_dic
              
#
#Read decays info in ENSDFZipFile
#            


## Read alpha decay info in ENSDF database for a given nucleus
# 
#  @param  Z atomic/charge number
#  @param A mass number
#     
def ReadAlphaDecayInENSDFZipFile(Z,A):
    
    daughter_Z=Z-2
    daughter_A=A-4
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    print "Look for alpha decays"
    print "daughter ",daughter_name
    
    qcalc_dic=QCalc.ReadQCALCResults(Z,A)
    parent_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A)
    daughter_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(daughter_Z,daughter_A)
          
    theZipFile=zipfile.ZipFile(zip_file_name_decays)
    list_names=theZipFile.namelist()
    if daughter_name in list_names:
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        in_a_block=False
        decay_list=[]
        ilevel=-1
        last_level_energy=0.
        for line in lines:
            print "Alpha ",line,line[5:9]
            if line[5:9].replace(" ","") == "":
                nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                    is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
                if (is_alpha_decay and not is_tentative):
                    in_a_block=True
                    decay_list+=[{}]
                    decay_list[-1]["dsid"]=dsid
                    decay_list[-1]["daughter_level_energy_list"]=[]
                    decay_list[-1]["alpha_energy_list"]=[]
                    decay_list[-1]["alpha_intensity_list"]=[]
                    decay_list[-1]["br"]=1.
                    ilevel=-1
                    if "Qalpha" in  qcalc_dic:
                        decay_list[-1]["Qval_qcalc"]=qcalc_dic["Qalpha"]
                    print "In a block",in_a_block
            if line=='\n':
                in_a_block=False
                #Check for zero intensity 
                if len(decay_list)>0:
                    energy_vec=numpy.array(decay_list[-1]["alpha_intensity_list"])
                    intensity_vec=numpy.array(decay_list[-1]["alpha_intensity_list"])
                    daughter_level_vec=numpy.array(decay_list[-1]["daughter_level_energy_list"])
                    indices=numpy.where(intensity_vec>0.)
                    if len(indices[0])>0:
                         decay_list[-1]["daughter_level_energy_list"]=daughter_level_vec[indices]
                         decay_list[-1]["alpha_energy_list"]=energy_vec[indices]
                         decay_list[-1]["alpha_intensity_list"]=intensity_vec[indices]
                    else:
                        if len(intensity_vec)>0:
                            decay_list[-1]["alpha_intensity_list"]=intensity_vec+100./len(intensity_vec)
                        else: #case where no decay info given altough alpha is given as potential decay, check with wallet cards
                            allowed_decays,decay_vec,Parent_elevel_list,Parent_half_life_list=NuclearWalletCards.ReadDecaysInWalletCards(Z,A)
                            print "test ",allowed_decays
                            if "Alpha" in allowed_decays:
                                print "How to do"
                            else:
                                decay_list=decay_list[0:-1]
                                
                                
                        
                        
                    
                    
                    
                    
            
            if in_a_block and line[5:9]=="  P " :
                nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given =ReadParentRecord(line)
                if (e_is_given):
                    index_plevel=FindLevel(parent_elevels_dic,e_level,ref_level,max_diff_e=2.)
                    if index_plevel >=0:
                        decay_list[-1]["parent_e_level"]=parent_elevels_dic["level_energies"][index_plevel]
                    else:
                        decay_list[-1]["parent_e_level"]=e_level
                else:
                    decay_list[-1]["parent_e_level"]=0.
                decay_list[-1]["parent_half_life"]=half_life
                decay_list[-1]["parent_half_life_str"]=half_life_str
                if is_Qval_given:
                    decay_list[-1]["Qval_ensdf"]=Qval
            if in_a_block and line[5:9]=="  N " :
                nr,nt,br,nb,np,is_nr_given,is_nt_given,is_br_given,is_nb_given,is_np_given=ReadNormalizationRecord(line)
                if is_br_given:
                    decay_list[-1]["br"]=br
            
            if in_a_block and line[5:9] =="  L "  and "FLAG" not in line:
                ilevel+=1
                energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,J_str=ReadLevelRecord(line)
                index_dlevel=FindLevel(daughter_elevels_dic,energy,ref_level,max_diff_e=2.)
                last_level_energy=energy
                if index_dlevel >=0:
                    last_level_energy=daughter_elevels_dic["level_energies"][index_dlevel]
            
            if in_a_block and line[5:9]=="  A " :
                energy,e_value_unit,intensity,hf,e_is_given,intensity_given,hf_given=ReadAlphaRecord(line)
                print "Alpha in ",line
                print energy,e_value_unit,intensity,hf,e_is_given,intensity_given,hf_given
                if e_is_given :
                    decay_list[-1]["alpha_energy_list"]+=[energy]
                    decay_list[-1]["alpha_intensity_list"]+=[intensity]
                    decay_list[-1]["daughter_level_energy_list"] +=[last_level_energy]  
                  
                
        if len(decay_list) >0:
            return  decay_list  
        else:
            print "No alpha decay for %s",daughter_name
            return None
    else:
            return None
## Read beta minus decay info in ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number
#  @param A mass number
#   
        
def ReadBminusDecayInENSDFZipFile(Z,A):
    
    daughter_Z=Z+1
    daughter_A=A
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    print "Look for beta minus"
    print "daugther ",daughter_name
    qcalc_dic=QCalc.ReadQCALCResults(Z,A)
   
    parent_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A)
    daughter_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(daughter_Z,daughter_A)

            
    theZipFile=zipfile.ZipFile(zip_file_name_decays)
    list_names=theZipFile.namelist()
    if daughter_name in list_names:
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        in_a_block=False
        parent_line_was_detected=False
        decay_list=[]
        ilevel=-1
        last_level_energy=0.
        for line in lines:
            if line[5:9].replace(" ","") == "":
                nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                    is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
                if (is_beta_minus_decay and not is_tentative):
                    in_a_block=True
                    decay_list+=[{}]
                    decay_list[-1]["dsid"]=dsid
                    decay_list[-1]["daughter_level_energy_list"]=[]
                    decay_list[-1]["betam_energy_list"]=[]
                    decay_list[-1]["betam_intensity_list"]=[]
                    decay_list[-1]["is_forbidden_decay_list"]=[]
                    decay_list[-1]["is_forbidden_unique_decay_list"]=[]
                    decay_list[-1]["forbidness_order_list"]=[]
                    decay_list[-1]["g4beta_type_flag_list"]=[]
                    decay_list[-1]["br"]=1.
                    ilevel=-1
                    if "Qbeta-" in  qcalc_dic:
                        decay_list[-1]["Qval_qcalc"]=qcalc_dic["Qbeta-"]
            if line=='\n':
                if (in_a_block and not parent_line_was_detected):
                    decay_list=decay_list[0:-1]
                in_a_block=False
                parent_line_was_detected=False
            
            if in_a_block and line[5:9]=="  P " :
                nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given =ReadParentRecord(line)
                parent_line_was_detected=True
                if (e_is_given):
                    index_plevel=FindLevel(parent_elevels_dic,e_level,ref_level,max_diff_e=2.)
                    if index_plevel >=0:
                        decay_list[-1]["parent_e_level"]=parent_elevels_dic["level_energies"][index_plevel]
                    else:
                        decay_list[-1]["parent_e_level"]=e_level
                else:
                    decay_list[-1]["parent_e_level"]=0.
                decay_list[-1]["parent_half_life"]=half_life
                decay_list[-1]["parent_half_life_str"]=half_life_str
                if is_Qval_given:
                    decay_list[-1]["Qval_ensdf"]=Qval
            if in_a_block and parent_line_was_detected and line[5:9]=="  N " :
                nr,nt,br,nb,np,is_nr_given,is_nt_given,is_br_given,is_nb_given,is_np_given=ReadNormalizationRecord(line)
                if is_br_given:
                    decay_list[-1]["br"]=br
                if is_nb_given:
                    decay_list[-1]["nb"]=nb
            
            if in_a_block and line[5:9] =="  L "  and "FLAG" not in line:
                ilevel+=1
                energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,J_str=ReadLevelRecord(line)
                index_dlevel=FindLevel(daughter_elevels_dic,energy,ref_level,max_diff_e=2.)
                last_level_energy=energy
                if index_dlevel >=0:
                    last_level_energy=daughter_elevels_dic["level_energies"][index_dlevel] 
            
            if in_a_block  and parent_line_was_detected and line[5:9]=="  B " :
                end_energy, e_value_unit,e_is_given,intensity,is_intensity_given,\
                   logft,is_logft_given,is_forbidden_decay,is_forbidden_unique_decay,forbidness_order, \
                      g4_beta_type_flag=ReadBetaMinusRecord(line)
                decay_list[-1]["betam_energy_list"]+=[end_energy]
                decay_list[-1]["betam_intensity_list"]+=[intensity]
                decay_list[-1]["daughter_level_energy_list"] +=[last_level_energy]
                decay_list[-1]["is_forbidden_decay_list"]+=[is_forbidden_decay]
                decay_list[-1]["is_forbidden_unique_decay_list"]+=[is_forbidden_unique_decay]
                decay_list[-1]["forbidness_order_list"]+=[forbidness_order]
                decay_list[-1]["g4beta_type_flag_list"]+=[g4_beta_type_flag]
                 
                
        if len(decay_list) >0:
            if len(decay_list[-1]["betam_energy_list"]) ==0:
                decay_list[-1]["betam_energy_list"]+=[-1.]
                decay_list[-1]["betam_intensity_list"]+=[decay_list[-1]["br"]]
                decay_list[-1]["daughter_level_energy_list"] +=[last_level_energy]
                decay_list[-1]["is_forbidden_decay_list"]+=[False]
                decay_list[-1]["is_forbidden_unique_decay_list"]+=[False]
                decay_list[-1]["forbidness_order_list"]+=[0]
                decay_list[-1]["g4beta_type_flag_list"]+=["allowed"]
            return  decay_list  
        else:
            print "No beta minus decay"
            return None
    else:
            print "No beta minus decay",daughter_name
            return None

## Read electron captutre  decay info in ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number 
#  @param A mass number
#    
def ReadECBplusDecayInENSDFZipFile(Z,A,is_xundl=False):
    
    #Get the daughter because radioactive ENSDF are specified in function of daughter
    daughter_Z=Z-1
    daughter_A=A
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    print "Look for beta plus"
    print "daugther ",daughter_name
    qcalc_dic=QCalc.ReadQCALCResults(Z,A)

    parent_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A)
    daughter_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(daughter_Z,daughter_A)

    
        
    theZipFile=zipfile.ZipFile(zip_file_name_decays)
    list_names=theZipFile.namelist()
    if daughter_name in list_names: #Check if the daughter has a file in the ensdf files
        #Read the lines of the ensdf file
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        if daughter_name in ['82RB','82KR'] and not is_xundl:
            file_object= theZipFile.open(daughter_name+"1")
            lines=file_object.readlines()
            file_object.close()

            
        in_a_block=False
        parent_line_was_detected=False
        decay_list=[]
        ilevel=-1
        last_level_energy=0.
        secondary_line_found=False
        for line in lines:
            if line[5:9].replace(" ","") == "":
                nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                    is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
                if (is_EC_beta_plus_decay and not is_tentative):
                    in_a_block=True
                    decay_list+=[{}]
                    decay_list[-1]["dsid"]=dsid
                    decay_list[-1]["daughter_level_energy_list"]=[]
                    decay_list[-1]["beta+_intensity_list"]=[]
                    decay_list[-1]["ec_intensity_list"]=[]
                    decay_list[-1]["ec_ck_intensity_list"]=[]
                    decay_list[-1]["ec_cm+_intensity_list"]=[]
                    decay_list[-1]["ec_cl_intensity_list"]=[]
                    decay_list[-1]["is_forbidden_decay_list"]=[]
                    decay_list[-1]["is_forbidden_unique_decay_list"]=[]
                    decay_list[-1]["forbidness_order_list"]=[]
                    decay_list[-1]["g4beta_type_flag_list"]=[]
                    decay_list[-1]["br"]=1.
                    ilevel=-1
                    if "QEC" in  qcalc_dic:
                        decay_list[-1]["Qval_EC_qcalc"]=qcalc_dic["QEC"]
                    if "Qbeta+" in  qcalc_dic:
                        decay_list[-1]["Qval_beta+_qcalc"]=qcalc_dic["Qbeta+"]
            
            #Check at the end of a block if a parent line was detected, if not it is rejected
            if line=='\n': 
                if (in_a_block and not parent_line_was_detected):
                    decay_list=decay_list[0:-1]
                in_a_block=False
                parent_line_was_detected=False
            
            if in_a_block and line[5:9]=="  P " :
                nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given =ReadParentRecord(line)
                parent_line_was_detected=True
                if (e_is_given):
                    index_plevel=FindLevel(parent_elevels_dic,e_level,ref_level,max_diff_e=2.)
                    if index_plevel >=0:
                        decay_list[-1]["parent_e_level"]=parent_elevels_dic["level_energies"][index_plevel]
                    else:
                        decay_list[-1]["parent_e_level"]=e_level
                else:
                    decay_list[-1]["parent_e_level"]=0.
                decay_list[-1]["parent_half_life"]=half_life
                decay_list[-1]["parent_half_life_str"]=half_life_str
                if is_Qval_given:
                    decay_list[-1]["Qval_ensdf"]=Qval
            if in_a_block and line[5:9]=="  N " and parent_line_was_detected:
                nr,nt,br,nb,np,is_nr_given,is_nt_given,is_br_given,is_nb_given,is_np_given=ReadNormalizationRecord(line)
                if is_br_given:
                    decay_list[-1]["br"]=br
                if is_nb_given:
                    decay_list[-1]["nb"]=nb
            
            if in_a_block and line[5:9] =="  L "  and "FLAG" not in line:
                ilevel+=1
                energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,J_str=ReadLevelRecord(line)
                index_dlevel=FindLevel(daughter_elevels_dic,energy,ref_level,max_diff_e=2.)
                last_level_energy=energy
                if index_dlevel >=0:
                    last_level_energy=daughter_elevels_dic["level_energies"][index_dlevel]    
            
            if in_a_block and parent_line_was_detected and line[5:9]=="  E " :
                energy, e_value_unit,e_is_given,IB,is_IB_given, \
                    IE,is_IE_given,TI,is_TI_given,logft,is_logft_given, \
                        is_forbidden_decay,is_forbidden_unique_decay,forbidness_order, \
                          g4_beta_type_flag=ReadECRecord(line)
                secondary_line_found=True
                print "Yes"
                    
                
                decay_list[-1]["beta+_intensity_list"]+=[IB]
                decay_list[-1]["ec_intensity_list"]+=[IE]
                decay_list[-1]["daughter_level_energy_list"] +=[last_level_energy] 
                decay_list[-1]["is_forbidden_decay_list"]+=[is_forbidden_decay]
                decay_list[-1]["is_forbidden_unique_decay_list"]+=[is_forbidden_unique_decay]
                decay_list[-1]["forbidness_order_list"]+=[forbidness_order]
                decay_list[-1]["ec_ck_intensity_list"]+=[0.]
                decay_list[-1]["ec_cm+_intensity_list"]+=[0.]
                decay_list[-1]["ec_cl_intensity_list"]+=[0.]
                decay_list[-1]["g4beta_type_flag_list"]+=[g4_beta_type_flag]
            if in_a_block and parent_line_was_detected and line[6:9]==" E "  and line[5] != " ":
                card_dic=ReadContinuationCard(line)
                if ("CK" in card_dic and'value' in  card_dic["CK"] ):
                    decay_list[-1]["ec_ck_intensity_list"][-1]=card_dic["CK"]['value'] 
                if ("CL" in card_dic and'value' in  card_dic["CL"] ):
                    decay_list[-1]["ec_cl_intensity_list"][-1]=card_dic["CL"]['value'] 
                if ("CM+" in card_dic and'value' in  card_dic["CM+"]):
                    decay_list[-1]["ec_cm+_intensity_list"][-1]=card_dic["CM+"]['value']  
                
                
        if len(decay_list) >0 and secondary_line_found:
            return  decay_list  
        else:
            print "No EC,beta+ decay"
            return None
    else:
            print "No EC beta+ decay"
            return None
        
## Check if an IT decay is fullx defined the ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number 
#  @param A mass number
#   
def CheckIfITDecayInENSDFZipFile(Z,A):
    daughter_Z=Z
    daughter_A=A
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    theZipFile=zipfile.ZipFile(zip_file_name_decays)
    list_names=theZipFile.namelist()
    if daughter_name in list_names:
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        in_a_block=False
        parent_line_was_detected=False
        decay_list=[]
        ilevel=-1
        last_level_energy=0.
        for line in lines:
            if line[5:9].replace(" ","") == "":
                nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                    is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
              
                if (is_IT_decay and not is_tentative):
                    in_a_block=True
                    decay_list+=[{}]
                    decay_list[-1]["dsid"]=dsid
                    decay_list[-1]["br"]=1.
                    ilevel=-1
            if line=='\n':
                if (in_a_block and not parent_line_was_detected):
                    decay_list=decay_list[0:-1]
                in_a_block=False
                parent_line_was_detected=False
                
            
            if in_a_block and line[5:9]=="  P " :
                nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given =ReadParentRecord(line)
                Jpi_vec=GetValuesForJPiField(Jpi_str)
                if len(Jpi_vec)==0:
                    Jpi_str="NA"
                return True,e_level,ref_level,half_life,Jpi_str
    return False,0.,0.,0.,""
    
## Read IT decay info in ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number 
#  @param A mass number
#       
def ReadITDecayInENSDFZipFile(Z,A):
    global adopted_levels_and_gammas,ids_isomer_to_add_in_database,elevels_isomer_to_add_in_database,elevels_daughter_isomer_to_add_in_database,half_life_isomer_to_add_in_database 

    daughter_Z=Z
    daughter_A=A
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    """
    print "Look for IT decay"
    print "daugther ",daughter_name
    """
    qcalc_dic=QCalc.ReadQCALCResults(Z,A)
  
    parent_elevels_dic=ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A)
   
    T_limit=1.e-9
   
    id_nuc=1000*A+Z
    iso_decay_list=[]
    decay_list=[]
    """
    search_array=numpy.array(ids_isomer_to_add_in_database)
    if id_nuc in ids_isomer_to_add_in_database:
        indices=numpy.where(search_array==id_nuc)
        for ind in indices[0]:
            decay_list+=[{}]
            decay_list[-1]["br"]=1.
            decay_list[-1]["parent_e_level"]=elevels_isomer_to_add_in_database[ind]
            half_life=half_life_isomer_to_add_in_database[ind]
            decay_list[-1]["parent_half_life"]=half_life
            decay_list[-1]["parent_half_life_str"]="%.4e" %(half_life)
    """
    it_energies=[]
        
        
        
    theZipFile=zipfile.ZipFile(zip_file_name_decays)
    list_names=theZipFile.namelist()
    if daughter_name in list_names:
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        in_a_block=False
        parent_line_was_detected=False
        ilevel=-1
        last_level_energy=0.
        for line in lines:
            if line[5:9].replace(" ","") == "":
                nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                    is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
              
                if (is_IT_decay and not is_tentative):
                    in_a_block=True
                    decay_list+=[{}]
                    decay_list[-1]["dsid"]=dsid
                    decay_list[-1]["br"]=1.
                    ilevel=-1
            if line=='\n':
                if (in_a_block and not parent_line_was_detected):
                    decay_list=decay_list[0:-1]
                in_a_block=False
                parent_line_was_detected=False
                
            
            if in_a_block and line[5:9]=="  P " :
                parent_line_was_detected=True
                nuc_name,e_level,e_is_given,ref_level,Jpi_str,half_life,half_life_str,is_stable,Qval,is_Qval_given =ReadParentRecord(line)
                if (e_is_given):
                    index_plevel=FindLevel(parent_elevels_dic,e_level,ref_level,max_diff_e=2.)
                    if index_plevel >=0:
                        decay_list[-1]["parent_e_level"]=parent_elevels_dic["level_energies"][index_plevel]
                    else:
                        decay_list[-1]["parent_e_level"]=e_level
                else:
                    decay_list[-1]["parent_e_level"]=0.     
                it_energies+=[decay_list[-1]["parent_e_level"]]
                
                
                decay_list[-1]["parent_half_life"]=half_life
                decay_list[-1]["parent_half_life_str"]=half_life_str
            
            if in_a_block and line[5:9]=="  N " :
                nr,nt,br,nb,np,is_nr_given,is_nt_given,is_br_given,is_nb_given,is_np_given=ReadNormalizationRecord(line)
                if is_br_given:
                    decay_list[-1]["br"]=br
        
    
    iso_decay_list=[]
    print parent_elevels_dic
    for i in range(len(parent_elevels_dic["level_energies"])):
        e=parent_elevels_dic["level_energies"][i]
        half_life=parent_elevels_dic["half_life_list"][i]
        if e >0 and e not in it_energies and half_life>T_limit:
            iso_decay_list+=[{}]
            iso_decay_list[-1]["br"]=1.
            iso_decay_list[-1]["parent_e_level"]=e

            iso_decay_list[-1]["parent_half_life"]=half_life
            iso_decay_list[-1]["parent_half_life_str"]="%.4e" %(half_life)
    if len(iso_decay_list)==0:
         iso_decay_list=None
    if len(decay_list)==0:
         decay_list=None
    print "iso_decay_list",iso_decay_list
    return decay_list,iso_decay_list    

## Read all decay infos in ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number
#  @param A mass number
#   
 

def ReadDecaysInENSDFZipFile(Z,A):
    global zip_file_name_decays, zip_file_name_levels, zip_file_name_decays_xundl

    nuc_name=NuclearWalletCards.GetNucleusName(Z, A)
    Parent_elevel_str_list=[]
    Parent_half_life_str_list=[]
    Parent_half_life_list=[]
    Parent_elevel_list=[]
    
    decays_dic={}
    allowed_decays=[]
    
    #Read the different decays
    ##########################
    
    BplusEC_decay_list=ReadECBplusDecayInENSDFZipFile(Z,A)
    IT_decay_list,ISO_decay_list=ReadITDecayInENSDFZipFile(Z,A,)
    Bminus_decay_list=ReadBminusDecayInENSDFZipFile(Z,A)
    Alpha_decay_list=ReadAlphaDecayInENSDFZipFile(Z,A)
    
    print "XUNDL"
    copy_zip_file_name_decays=zip_file_name_decays
    zip_file_name_decays=zip_file_name_decays_xundl
    BplusEC_decay_list_xundl=ReadECBplusDecayInENSDFZipFile(Z,A,is_xundl=True)
    IT_decay_list_xundl,ISO_decay_list_xundl=ReadITDecayInENSDFZipFile(Z,A)
    Bminus_decay_list_xundl=ReadBminusDecayInENSDFZipFile(Z,A,)
    Alpha_decay_list_xundl=ReadAlphaDecayInENSDFZipFile(Z,A)
    zip_file_name_decays=copy_zip_file_name_decays
    
    
    
    
    decays_vec=[BplusEC_decay_list,Bminus_decay_list,Alpha_decay_list,IT_decay_list,ISO_decay_list]
    decays_vec_xundl=[BplusEC_decay_list_xundl,Bminus_decay_list_xundl,
                      Alpha_decay_list_xundl,IT_decay_list_xundl,ISO_decay_list_xundl]
    name_decays=["BetaPlusEC","BetaMinus","Alpha","IT","ISO"]
    
    for i in range(len(decays_vec)):
        if (decays_vec[i] is  None or len(decays_vec[i])<=0):
            decays_vec[i]=decays_vec_xundl[i]
        
    for i in range(len(decays_vec)):
        if (decays_vec[i] is not None and len(decays_vec[i])>0):
            map_pelevel_str_to_decay_index={}
            elevel_vec,new_decay_list,parent_half_life_vec,parent_half_life_str_vec=OrderDecaysWithIncreasingParentEnergy(decays_vec[i])
            decays_dic[name_decays[i]]={"decay_list":new_decay_list}
            allowed_decays+=[name_decays[i]]
            j=-1
            for e in elevel_vec:
                j+=1
                estr="%.6f" %(e)
                parent_half_life_str=parent_half_life_str_vec[j]
                if estr not in  Parent_elevel_str_list:
                    Parent_elevel_str_list+=[estr]
                    Parent_elevel_list+=[e]
                    Parent_half_life_str_list+=[parent_half_life_str]
                    Parent_half_life_list+=[parent_half_life_vec[j]]
                if estr not in map_pelevel_str_to_decay_index:
                    map_pelevel_str_to_decay_index[estr]=[]
                map_pelevel_str_to_decay_index[estr]+=[j]
            decays_dic[name_decays[i]]["map_pelevel_str_to_decay_index"]=map_pelevel_str_to_decay_index
    
    if (len(allowed_decays)==0 or allowed_decays==["ISO"]):
        decays_vec=[BplusEC_decay_list_xundl,Bminus_decay_list_xundl,Alpha_decay_list_xundl,IT_decay_list_xundl,ISO_decay_list_xundl]
        name_decays=["BetaPlusEC","BetaMinus","Alpha","IT","ISO"]
        decays_dic={}
        allowed_decays=[]
       
        for i in range(len(decays_vec)):
            if (decays_vec[i] is not None and len(decays_vec[i])>0):
                map_pelevel_str_to_decay_index={}
                elevel_vec,new_decay_list,parent_half_life_vec,parent_half_life_str_vec=OrderDecaysWithIncreasingParentEnergy(decays_vec[i])
                decays_dic[name_decays[i]]={"decay_list":new_decay_list}
                allowed_decays+=[name_decays[i]]
                j=-1
                for e in elevel_vec:
                    j+=1
                    estr="%.6f" %(e)
                    parent_half_life_str=parent_half_life_str_vec[j]
                    if estr not in  Parent_elevel_str_list:
                        Parent_elevel_str_list+=[estr]
                        Parent_elevel_list+=[e]
                        Parent_half_life_str_list+=[parent_half_life_str]
                        Parent_half_life_list+=[parent_half_life_vec[j]]
                    if estr not in map_pelevel_str_to_decay_index:
                        map_pelevel_str_to_decay_index[estr]=[]
                    map_pelevel_str_to_decay_index[estr]+=[j]
                    decays_dic[name_decays[i]]["map_pelevel_str_to_decay_index"]=map_pelevel_str_to_decay_index
   
    
    if (len(allowed_decays)>0):
        Parent_elevel_str_list=numpy.array(Parent_elevel_str_list)
        Parent_elevel_list=numpy.array(Parent_elevel_list)
        Parent_half_life_str_list=numpy.array(Parent_half_life_str_list)
        Parent_half_life_list=numpy.array(Parent_half_life_list)
        indices=numpy.argsort(Parent_elevel_list)
        Parent_elevel_str_list=Parent_elevel_str_list[indices]
        Parent_elevel_list=Parent_elevel_list[indices]
        Parent_half_life_str_list=Parent_half_life_str_list[indices]
        Parent_half_life_list=Parent_half_life_list[indices]        
    
     
    return nuc_name, allowed_decays, decays_dic, Parent_elevel_str_list, Parent_elevel_list, Parent_half_life_str_list, Parent_half_life_list

## Write list of IT decays found in  ENSDF database 
# 
#
def WriteListOfITDecays(file_name="../isomers_ENDSF_August2012/list_Isomers_in_ENSDF_HalfLife_gt_1e-9_new.txt",
                        file_ground_state_name="../isomers_ENDSF_August2012/list_ground_states_in_ENSDF.txt"):
    Zmin=1
    Zmax=120
    decay_list,level_list=ReadListOfAvailableNucleiInENDSF()
    map_AZ_to_names, map_names_to_AZ,A_vec,Z_vec = NuclearWalletCards.ReadNucleusNameList()
   
    header_line="Name\t\tZ\t\tA\t\tElevel[keV]\t\tT1/2[s]\t\t\t\tJpi\t\tMag Dipole Moment[magneton]\t\tENSDF type\n"
    file_obj=open(file_name,'w')
    file_obj.write(header_line)
    file_obj.close()
    lines_f_object=""
    file_obj2=None
    last_Z=-1
    file_obj3=open(file_ground_state_name,'w')
    file_obj3.write(header_line)
    lines_ground_state=""
    file_obj3.close()
    for i in range(len(Z_vec)):
        Z=Z_vec[i]
        A=A_vec[i]
        nuc_name=NuclearWalletCards.GetNucleusName(Z,A)
        print Z,A
        if Z>=Zmin and Z<=Zmax:
            if Z>last_Z:
                if file_obj2 is not None:
                    file_obj2.close()
                file_obj=open(file_name,'a')
                file_obj.write(lines_f_object)
                file_obj.close()
                lines_f_object=""
                file_obj3=open(file_ground_state_name,'a')
                file_obj3.write(lines_ground_state)
                file_obj3.close()
                
                lines_ground_state=""
                file_obj2=open("../isomers_ENDSF_August2012/list_Isomers_in_ENSDF_Z%i_new.txt" %(Z),'w')
                file_obj2.write(header_line)
            last_Z=Z
            res,elevel,ref_level,half_life,Jpi_str=CheckIfITDecayInENSDFZipFile(Z,A,zip_file_name="../ensdf_data/august_2012/all_decays.zip")
            if res :
                str="%6s\t\t%3i\t\t%3i\t\t%.4e\t\t%.4e\t%12s\t%12s" %(nuc_name,Z,A,elevel,half_life,Jpi_str,"NA") 
                if ref_level!="Zero":
                    str="%6s\t\t%3i\t\t%3i\t\t%.4e+%s\t\t%.4e\t%12s\t%12s" %(nuc_name,Z,A,elevel,ref_level,half_life,Jpi_str,"NA")
                str+="\t\t\t\t\t\tIT\n"
                file_obj2.write(str)
                lines_f_object+=str
            level_energy_vec,level_half_life_vec,level_is_stable_vec,ref_level_vec,level_magnetic_dipole_moment_str,level_Jpi_str=ReadAdoptedLevelInENSDFZipFile(Z,A,zip_file_name="../ensdf_data/august_2012/all_levels.zip")
            if level_energy_vec is not None:
                for i in range(len(level_energy_vec)):
                    elevel=level_energy_vec[i]
                    ref_level=ref_level_vec[i]
                    half_life=level_half_life_vec[i]
                    mmom_str=level_magnetic_dipole_moment_str[i].replace(" ","")
                    Jpi_str=level_Jpi_str[i].replace(" ","")
                    if (elevel>0 or ref_level!="Zero"):
                        str="%6s\t\t%3i\t\t%3i\t\t%.4e\t\t%.4e\t%12s\t%12s" %(nuc_name,Z,A,elevel,half_life,Jpi_str,mmom_str) 
                        if ref_level!="Zero":
                            str="%6s\t\t%3i\t\t%3i\t\t%.4e+%s\t\t%.4e\t%12s\t%12s" %(nuc_name,Z,A,elevel,ref_level,half_life,Jpi_str,mmom_str)
                        str+="\t\t\t\t\t\tALG\n"
                        file_obj2.write(str)
                        if (half_life>=1.e-9):
                            lines_f_object+=str
                    
                    if (elevel==0. and ref_level=="Zero"):
                        str="%6s\t\t%3i\t\t%3i\t\t%.4e\t\t%.4e\t%12s\t%12s" %(nuc_name,Z,A,elevel,half_life,Jpi_str,mmom_str) 
                        if level_is_stable_vec[i]:
                            str="%6s\t\t%3i\t\t%3i\t\t%.4e\t\t%10s\t%12s\t%12s" %(nuc_name,Z,A,elevel,"STABLE",Jpi_str,mmom_str) 
                        str+="\t\t\t\t\t\tALG\n"
                        lines_ground_state+=str
                        
                        
            
    if file_obj2 is not None:
        file_obj2.close()
    file_obj=open(file_name,'a')
    file_obj.write(lines_f_object)
    file_obj.close()
    
    file_obj3=open(file_ground_state_name,'a')
    file_obj3.write(lines_ground_state)
    file_obj3.close()

## Read adopted level info in ENSDF database for a given nucleus
# 
#  @param  Z atomic/charge number
#  @param A mass number
#  
def ReadAdoptedLevelInENSDFZipFile(Z,A):
    daughter_Z=Z
    daughter_A=A
    daughter_name=NuclearWalletCards.GetNucleusName(daughter_Z, daughter_A)
    theZipFile=zipfile.ZipFile(zip_file_name_levels)
    list_names=theZipFile.namelist()
    if daughter_name in list_names:
        file_object= theZipFile.open(daughter_name)
        lines=file_object.readlines()
        file_object.close()
        level_energy_vec=[]
        level_half_life_vec=[]
        level_is_stable_vec=[]
        level_magnetic_dipole_moment_str=[]
        level_Jpi_str=[]
        
        
        ref_level_vec=[]
        if lines[0][5:9].replace(" ","") == "":
            nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(lines[0])
            
            if (is_adopted_levels and not is_tentative and not is_not_observed):
                for line in lines: 
                    if line[5:9] =="  L ":
                        level_energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,Jpi_str=ReadLevelRecord(line)
                        level_energy_vec+=[level_energy]
                        level_half_life_vec+=[half_life]
                        ref_level_vec+=[ref_level]
                        level_magnetic_dipole_moment_str+=["NA"]
                        level_Jpi_str+=["NA"]
                        level_is_stable_vec+=[is_stable]
                        if (len(JPi_vec)) >0:
                            level_Jpi_str[-1]=Jpi_str
                            
                            
                        
                    if line[6:9] ==" L "  and line[5] !=" ":
                        card_dic=ReadContinuationCard(line)
                        if "MOMM1" in card_dic:
                            if "value_str" in card_dic["MOMM1"]:
                                level_magnetic_dipole_moment_str[-1]=card_dic["MOMM1"]["value_str"]
                                
                            
                            
                    
        return numpy.array(level_energy_vec),numpy.array(level_half_life_vec),level_is_stable_vec,ref_level_vec,level_magnetic_dipole_moment_str,level_Jpi_str
    else:
        return None,None,None,None,None,None
    
def FindLevel(levels_dic,level_energy,level_ref,max_diff_e=2.):
    index_level=-1
    if "level_reference_list" in levels_dic and len(levels_dic["level_reference_list"])>0. :
        indices=numpy.where(numpy.array(levels_dic["level_reference_list"]) == level_ref)
        if len(indices[0])>0:
            diff_array=numpy.abs(numpy.array(levels_dic["ensdf_level_energies"])[indices]-level_energy)
            arg_min=numpy.argmin(diff_array)
            if (diff_array[arg_min]<max_diff_e):
                index_level=indices[0][arg_min]
    return index_level

    
    
    
## Read adopted level and gamma infos in ENSDF database for a given nucleus
# 
#  @param Z atomic/charge number
#  @param A mass number
#   
def ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A,add_info_from_decay_file=True):
    de_pro_ref_level=5.
    nb_ref_levels=1
    
    nuc_name=NuclearWalletCards.GetNucleusName(Z, A)
    if nuc_name in adopted_levels_and_gammas:
        return adopted_levels_and_gammas[nuc_name]
    
    #Get level and gammas from the level file
    de_pro_ref_level=5.
    nb_ref_levels=1
    adopted_level_and_gammas_dic={"ensdf_level_energies":[],
                                  "level_reference_list":[],
                                  "half_life_list":[],"gamma_energies":[],
                                  "gammas":[],
                                  "JPi_vec_list":[],
                                  "ReferenceLevels":{"Zero":0,"U":1.e-6,"V":2.e-6,"W":3.e-6,
                                                     "X":4.e-6,"Y":5.e-6,"Z":6.e-6}}
    ReferenceLevelsAlsoInDecayFile=["Zero"]
  
    theZipFile=zipfile.ZipFile(zip_file_name_levels)
    list_names=theZipFile.namelist()
    if nuc_name in list_names:
        file_object= theZipFile.open(nuc_name)
        lines=file_object.readlines()
        file_object.close()
        valid_gamma=False
        g_energies=[]
        g_intensities=[]
        cc_vec=[]
        last_level_energy=-1.
        if lines[0][5:9].replace(" ","") == "":
            nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(lines[0])
            
            if (is_adopted_levels and not is_tentative and not is_not_observed):
                for line in lines: 
                    print line
                    if line[5:9] =="  L ":
                        valid_gamma=False
                        level_energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,J_str=ReadLevelRecord(line)
                        #print level_energy,e_is_given
                        if (e_is_given):
                            last_level_energy=level_energy
                            adopted_level_and_gammas_dic["ensdf_level_energies"]+=[level_energy]
                            adopted_level_and_gammas_dic["level_reference_list"]+=[ref_level]
                            adopted_level_and_gammas_dic["half_life_list"]+=[half_life]
                            adopted_level_and_gammas_dic["JPi_vec_list"]+=[JPi_vec]
                            adopted_level_and_gammas_dic["gamma_energies"]+=[[]]
                            adopted_level_and_gammas_dic["gammas"]+=[[]]
                            #This is added for +X,Y,Z levels
                            if ref_level not in adopted_level_and_gammas_dic["ReferenceLevels"]:
                                adopted_level_and_gammas_dic["ReferenceLevels"][ref_level]=nb_ref_levels*de_pro_ref_level
                                nb_ref_levels+=1
                                                  
                                
                        g_energies=[]
                        g_intensities=[]
                        cc_vec=[]
               
                        
                    if line[5:9] =="  G ":
                        valid_gamma=False
                        energy, e_value_unit,e_is_given,ref_level,RI,is_RI_given, \
                                MP_str,MR,is_MR_given,CC, is_CC_given,TI,is_TI_given \
                                                                    =ReadGammaRecord(line)
                                                                    
                        if (e_is_given and len(adopted_level_and_gammas_dic["ensdf_level_energies"])>1):
                            
                           
                            adopted_level_and_gammas_dic["gamma_energies"][-1]+=[energy]
                            gamma={"card":line,"ref_level":ref_level}
                            g_energies+=[energy]
                            valid_gamma=True
                            if len(MP_str)>0:
                                gamma["MP"]=MP_str
                            if (is_MR_given):
                                gamma["MR"]=MR
                            if (is_RI_given):
                                gamma["RI"]=RI
                                g_intensities+=[RI]
                            elif (is_CC_given and  is_TI_given):
                                RI=TI/(1.+CC)
                                gamma["RI"]=RI
                                g_intensities+=[RI]
                            else:
                                g_intensities+=[0.00000000]
                                
                            if (is_CC_given):
                                gamma["CC"]=CC
                                cc_vec+=[CC]
                            else:
                                cc_vec+=[0.0]
                            if (is_TI_given):
                                gamma["TI"]=TI
                            adopted_level_and_gammas_dic["gammas"][-1]+=[gamma]
                    if line[6:9] ==" G "  and line[5] !=" " and valid_gamma:
                        card_dic=ReadContinuationCard(line)
                        """
                        if "FL" in card_dic:
                            print "FL",card_dic["FL"],adopted_level_and_gammas_dic["ensdf_level_energies"][-1],adopted_level_and_gammas_dic["level_reference_list"]
                            print ""
                            print ""
                            print "\n\n\n"
                        """
                        for key in card_dic:
                            if key == "CC":
                                gamma=adopted_level_and_gammas_dic["gammas"][-1][-1]
                                if "CC" not in gamma and "value" in card_dic:
                                    gamma["CC"]=card_dic["value"]
                            elif key == "RI":
                                gamma=adopted_level_and_gammas_dic["gammas"][-1][-1]
                                if "RI" not in gamma and "value" in card_dic:
                                    gamma["RI"]=card_dic["value"]
                            else:
                                adopted_level_and_gammas_dic["gammas"][-1][-1][key]=card_dic[key]
                                
                                                                   
    list_of_buggy_nuclei=["156PM"]

    if (add_info_from_decay_file and nuc_name not in list_of_buggy_nuclei):
        theZipFile=zipfile.ZipFile(zip_file_name_decays)
        list_names=theZipFile.namelist()
        if nuc_name in list_names:
            file_object= theZipFile.open(nuc_name)
            lines=file_object.readlines()
            file_object.close()
            is_tentative=False
            is_not_observed=False
            index_elevel=-1
            index_gamma=-1
            g_energies=[]
            g_intensities=[]
            L_already_detected=False
            valid_gamma=False
            for line in lines:
                if line[5:9].replace(" ","") == "":
                    nuc_id,dsid, is_beta_minus_decay, is_EC_beta_plus_decay,is_alpha_decay,\
                        is_IT_decay,is_adopted_levels,is_gammas,is_tentative, is_not_observed\
                                                                = ReadIdentificationRecord(line)
                    L_already_detected=False
                    g_energies=[]
                    g_intensities=[]
                    index_elevel=-1
                    index_gamma=-1
                    g_energies=[]
                    g_intensities=[]
                    L_already_detected=False
                    
                if line[5:9] =="  L " and not is_tentative and not is_not_observed:
                    #
                    #fill the precedent level with gamma RI if not alraedy in
                    #
                    if len(g_intensities)>0 and L_already_detected:
                        g_intensities=numpy.array(g_intensities)
                        max_intensity=numpy.max(g_intensities)
                        if max_intensity>0.:
                            g_intensities=100.*g_intensities/max_intensity
                        for j in range(len(g_energies)):
                            e=g_energies[j]
                            intensity=g_intensities[j]
                            diff_array=numpy.abs(numpy.array(adopted_level_and_gammas_dic["gamma_energies"][index_elevel])-e)
                            index_gamma=numpy.argmin(diff_array)
                            gamma=adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]
                            if ("RI" not in gamma or gamma["RI"]<=0):
                                gamma["RI"]=intensity
                                adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]=gamma
                    
                    g_energies=[]
                    g_intensities=[]
                    
                    level_energy,e_value_unit,half_life,e_is_given,ref_level,is_stable,JPi_vec,meta_str,J_str=ReadLevelRecord(line)
                    if ref_level not in ReferenceLevelsAlsoInDecayFile:
                            ReferenceLevelsAlsoInDecayFile+=[ref_level]
                     
                    if (e_is_given):
                        L_already_detected=True
                        #Check if level does not exist already
                        new_level=False
                        index_elevel=FindLevel(adopted_level_and_gammas_dic,level_energy,ref_level,max_diff_e=2.)
                        if (index_elevel<0):
                            L_already_detected=False
                            """
                            new_level=True
                            adopted_level_and_gammas_dic["ensdf_level_energies"]+=[level_energy]
                            adopted_level_and_gammas_dic["half_life_list"]+=[half_life]
                            adopted_level_and_gammas_dic["JPi_vec_list"]+=[JPi_vec]
                            adopted_level_and_gammas_dic["gamma_energies"]+=[[]]
                            adopted_level_and_gammas_dic["gammas"]+=[[]]
                            adopted_level_and_gammas_dic["level_reference_list"]+=[ref_level]
                            if ref_level not in adopted_level_and_gammas_dic["ReferenceLevels"]:
                                adopted_level_and_gammas_dic["ReferenceLevels"][ref_level]=nb_ref_levels*de_pro_ref_level
                                nb_ref_levels+=1
                            index_elevel=-1
                            """
                            
                    else:
                        L_already_detected=False   
                            
                if line[5:9] =="  G " and L_already_detected:
                    valid_gamma=False
                    energy, e_value_unit,e_is_given,ref_level,RI,is_RI_given, \
                            MP_str,MR,is_MR_given,CC, is_CC_given,TI,is_TI_given \
                                                                    =ReadGammaRecord(line)
                    if  (e_is_given):
                        new_gamma=False
                        if len(adopted_level_and_gammas_dic["gamma_energies"])>0 and len(adopted_level_and_gammas_dic["gamma_energies"][index_elevel])>0:
                            diff_array=numpy.abs(numpy.array(adopted_level_and_gammas_dic["gamma_energies"][index_elevel])-energy)
                            arg_min=numpy.argmin(diff_array)
                            if (diff_array[arg_min]<2.):
                                index_gamma=arg_min
                            else:
                                new_gamma=True
                        else:
                            new_gamma=True
                        gamma={"card":line}
                        if (new_gamma):
                            index_gamma=-1
                            adopted_level_and_gammas_dic["gamma_energies"][index_elevel]+=[energy]
                            adopted_level_and_gammas_dic["gammas"][index_elevel]+=[gamma]
                        
                        gamma=adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]
                        valid_gamma=True
                        g_energies+=[energy]
                        if (is_RI_given):
                            g_intensities+=[RI]
                        elif (is_CC_given and  is_TI_given):
                            RI=TI/(1.+CC)
                            g_intensities+=[RI]
                        else:
                            g_intensities+=[0.00000000]
                        
                        if (is_CC_given and "CC" not in gamma.keys()):
                            gamma["CC"]=CC
                        if (is_TI_given and "TI" not in gamma.keys()):
                            gamma["TI"]=TI
                            
                        adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]=gamma
                        
                        
                        
                if line[6:9] ==" G " and line[5]!=" " and  valid_gamma and L_already_detected:
                        card_dic=ReadContinuationCard(line)
                        for key in card_dic:
                            if key not in adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]:
                                adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma][key]=card_dic[key]
                                
            index_elevel=-1
            if len(g_intensities)>0:
                    g_intensities=numpy.array(g_intensities)
                    g_intensities=100.*g_intensities/numpy.max(g_intensities)
                    for j in range(len(g_energies)):
                        e=g_energies[j]
                        intensity=g_intensities[j]
                        diff_array=numpy.abs(numpy.array(adopted_level_and_gammas_dic["gamma_energies"][index_elevel])-e)
                        if len(diff_array)>0:
                            index_gamma=numpy.argmin(diff_array)
                            gamma=adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]
                            if ("RI" not in gamma or gamma["RI"]<=0):
                                gamma["RI"]=intensity
                                adopted_level_and_gammas_dic["gammas"][index_elevel][index_gamma]=gamma
            
    
    
    #
    #Correct gamma energy if needed from F level information
    #
    adopted_level_and_gammas_dic["level_energies"]=[]
    for i in range(len(adopted_level_and_gammas_dic["ensdf_level_energies"])):
        p_level_ekin=adopted_level_and_gammas_dic["ensdf_level_energies"][i]
        p_level_ref=adopted_level_and_gammas_dic["level_reference_list"][i]
        p_level_ekin+=adopted_level_and_gammas_dic["ReferenceLevels"][p_level_ref]
        """
        if p_level_ref not in ReferenceLevelsAlsoInDecayFile and p_level_ref!="Zero":
            p_level_ekin=-9999999.
        """
        adopted_level_and_gammas_dic["level_energies"]+=[p_level_ekin]
        for j in range(len(adopted_level_and_gammas_dic["gamma_energies"][i])):
            gamma_energy=adopted_level_and_gammas_dic["gamma_energies"][i][j]
            gamma=adopted_level_and_gammas_dic["gammas"][i][j]
            if ("FL" in gamma.keys()):
                val_str=gamma["FL"]["value_str"]
                value,is_given,value_type,value_unit,directly_measured,ref_level =GetValueForEField(val_str)
                if (is_given):
                    final_level_index=FindLevel(adopted_level_and_gammas_dic,value,ref_level,max_diff_e=2.)
                    if (final_level_index>=0):
                        f_level_ekin=adopted_level_and_gammas_dic["ensdf_level_energies"][final_level_index]
                        f_level_ref=adopted_level_and_gammas_dic["level_reference_list"][final_level_index]
                        f_level_ekin+=adopted_level_and_gammas_dic["ReferenceLevels"][f_level_ref]
                        recomputed_gamma_energy=p_level_ekin-f_level_ekin
                        if (recomputed_gamma_energy>0  and abs(gamma_energy-recomputed_gamma_energy)>2.):
                            adopted_level_and_gammas_dic["gamma_energies"][i][j]=recomputed_gamma_energy
                            
    #
    #Sort the levels
    #
    indices=numpy.argsort(adopted_level_and_gammas_dic["level_energies"])
    elevels=[]
    for index in indices:
        e_level=adopted_level_and_gammas_dic["level_energies"][index]
        p_level_ref=adopted_level_and_gammas_dic["level_reference_list"][index]
        """
        if p_level_ref not in ReferenceLevelsAlsoInDecayFile:
            e_level=-99999999.
        """
        elevels+=[e_level]
        
    for key in adopted_level_and_gammas_dic:
        if key not in ["ReferenceLevels"]:
            new_list=[]
            i=0
            for index in indices:
                e_level=elevels[i]
                if e_level>=0.:
                    new_list+=[adopted_level_and_gammas_dic[key][index]]
                i+=1
            adopted_level_and_gammas_dic[key]=new_list
    
  
    adopted_levels_and_gammas[nuc_name]=adopted_level_and_gammas_dic             
    return adopted_level_and_gammas_dic
                       
## Order decays following increasing parent energy
# 
#
def OrderDecaysWithIncreasingParentEnergy(decays):
    parent_elevel_vec=[]
    parent_half_life_vec=[]
    parent_half_life_str_vec=[]
    new_decays=decays
    if len(decays) >0:
        for i in range(len(decays)):
            parent_elevel_vec+=[decays[i]["parent_e_level"]]
            parent_half_life_vec+=[decays[i]["parent_half_life"]]
            parent_half_life_str_vec+=[decays[i]["parent_half_life_str"]]
        parent_elevel_vec=numpy.array(parent_elevel_vec)
        parent_half_time_vec=numpy.array(parent_half_life_vec)
        parent_half_time_str_vec=numpy.array(parent_half_life_str_vec)
        arg_sort=numpy.argsort(parent_elevel_vec)
        new_decays=numpy.array(decays)[arg_sort]
        parent_elevel_vec=parent_elevel_vec[arg_sort]
        parent_half_time_vec=parent_half_time_vec[arg_sort]
        parent_half_time_str_vec=parent_half_time_str_vec[arg_sort]
    return parent_elevel_vec,new_decays,parent_half_life_vec,parent_half_life_str_vec





       
    
    
        
    






    
    
    
    

     