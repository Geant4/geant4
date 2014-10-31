## @package G4Database
#  G4Database PYTHON module to produce G4 RadDEcay and PhotoEvaporation database
#
# 


#author:  L.Desorgher
#
#History: 
#---------
#      21/07/2014     Creation


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
import Bricc 
import  ENSDF 

uma_in_MeV=931.494028

#
#Functions to write the database
#

## Write the file of the G4RadioactiveDecay database for a given nucleus
#
#  @param Z atomic/charge number
#  @param A mass number
#  @param new_directory directory where the new file with name z*.a* will be set
#
def WriteG4RadDecayFile(Z,
                        A,new_directory="../RadioactiveDecay4.0"):
    
    buggy_nuclei=["162TM"]
    nuc_name=NuclearWalletCards.GetNucleusName(Z, A)
    if nuc_name in buggy_nuclei:
        return 
    
    pstring_template="P   388.532000   1.0090e+04"
    string_with_20_blanks="                    "
    string_with_40_blanks=string_with_20_blanks+string_with_20_blanks
    nuc_name, allowed_decays, decays_dic, Parent_elevel_str_list, Parent_elevel_list, \
                           Parent_half_life_str_list, Parent_half_life_list = ENSDF.ReadDecaysInENSDFZipFile(Z,A)
    
    
   
    
    lines=[]
    allowed_decays1,decay_vec1,Parent_elevel_list1,Parent_half_life_list1=NuclearWalletCards.ReadDecaysInWalletCards(Z,A)

    
    
    
    if len(allowed_decays)==0 or allowed_decays==["ISO"]:
        print "allowed_decays1",allowed_decays1
        if len(allowed_decays1)>0:
            lines="# %s ( %.4e s)\n" %(nuc_name,Parent_half_life_list1[0])
            lines+="#  Excitation   Halflife Mode       Daughter Ex  Intensity           Q\n"
            for i in range(len(decay_vec1)):
                elevel_str="%.6f" %(Parent_elevel_list1[i])
                parent_elevel=Parent_elevel_list1[i]
                half_life_str="%.4e" %(Parent_half_life_list1[i])
                lines+="P"+string_with_20_blanks[0:13-len(elevel_str)]+elevel_str
                lines+=string_with_20_blanks[0:13-len(half_life_str)]+half_life_str+"\n"
                decay=decay_vec1[i]
                k=0
                print decay
                for channel in decay['type_list']:
                    lines+=string_with_40_blanks[0:34-len(channel)]+channel+"       0.0000"+"   %.4e"%(decay['BR_list'][k]/100.)+"\n"
                    k+=1
                k=0
                for channel in decay['type_list']:
                    intensity=decay['BR_list'][k]
                    elevel_str="%.4f" %(0.)
                    elevel_str=string_with_20_blanks[0:13-len(elevel_str)]+elevel_str
                    intensity_str="%.4e" %(intensity)
                    intensity_str=string_with_20_blanks[0:13-len(intensity_str)]+intensity_str
                    qval_str="%.4f" %(decay['Q_list'][k])
                    qval_str=string_with_20_blanks[0:13-len(qval_str)]+qval_str
                    k+=1
                    if intensity>0 :
                        lines+=string_with_40_blanks[0:34-len(channel)]+channel+elevel_str+intensity_str+qval_str+"\n"
                
        
    if len(allowed_decays)>0 :
        lines="# %s ( %s)\n" %(nuc_name,Parent_half_life_str_list[0])
        lines+="#  Excitation   Halflife Mode       Daughter Ex  Intensity           Q\n"
        for i in range(len(Parent_elevel_str_list)):
            elevel_str=Parent_elevel_str_list[i]
            parent_elevel=Parent_elevel_list[i]
            half_life_str="%.4e" %(Parent_half_life_list[i])
            lines+="P"+string_with_20_blanks[0:13-len(elevel_str)]+elevel_str
            lines+=string_with_20_blanks[0:13-len(half_life_str)]+half_life_str+"\n"
            
            
            #Complete with Wallet cards info
            #Serach the level index in wallet card decay
            decay_wallet=None
            if len(Parent_elevel_list1):
                find_vec=np.absolute(np.array(Parent_elevel_list1)-parent_elevel)
                index_level_wallet=np.argmin(find_vec)
                if find_vec[index_level_wallet]<2.:
                    decay_wallet=decay_vec1[index_level_wallet]
            
            #check the type of decay channel available from this parent elevel and get the branching ratios
            decay_names=["Alpha","BetaMinus","BetaPlusEC","IT","ISO"]
            name_occuring_decays=[]
            name_occuring_decay_channels=[] #need to add channel of BetaPlusEC
            occuring_decays_for_channels=[]
            br_vec=[]
            br_channels_vec=[] #need to add channel of BetaPlusEC
            
            for name in decay_names:
                if name=="BetaPlusEC":
                    name_wallet="BetaPlus"
                if name in decays_dic:
                    decay_list=decays_dic[name]["decay_list"]
                    map_pelevel_str_to_decay_index=decays_dic[name]["map_pelevel_str_to_decay_index"]
                    if (elevel_str in map_pelevel_str_to_decay_index):
                        name_occuring_decays+=[name]
                        index=map_pelevel_str_to_decay_index[elevel_str][0]
                        decay=decay_list[index]
                        br=decay['br']
                        br_vec+=[br]
                        if name !="BetaPlusEC":
                            if name != "ISO" or len(name_occuring_decay_channels)==0:
                                name_occuring_decay_channels+=[name]
                                br_channels_vec+=[br]
                                occuring_decays_for_channels+=[decay]
                        else:
                            nb=1.
                            if "nb" in decay:
                                nb=decay["nb"]
                            ec_cm_rel_intensities=numpy.array(decay['ec_cm+_intensity_list'])
                            ec_cl_rel_intensities=numpy.array(decay['ec_cl_intensity_list'])
                            ec_ck_rel_intensities=numpy.array(decay['ec_ck_intensity_list'])
                            total_rel_intensities=ec_cm_rel_intensities+ec_cl_rel_intensities+ec_ck_rel_intensities
                            indices=numpy.where(total_rel_intensities>0)
                            ec_cm_rel_intensities[indices]=ec_cm_rel_intensities[indices]/total_rel_intensities[indices]
                            ec_cl_rel_intensities[indices]=ec_cl_rel_intensities[indices]/total_rel_intensities[indices]
                            ec_ck_rel_intensities[indices]=ec_ck_rel_intensities[indices]/total_rel_intensities[indices]
                            
                            ec_intensities=numpy.array(decay['ec_intensity_list'])
                            ec_intensity=sum(ec_intensities)
                            ec_m_intensities=ec_cm_rel_intensities*ec_intensities
                            ec_l_intensities=ec_cl_rel_intensities*ec_intensities
                            ec_k_intensities=ec_ck_rel_intensities*ec_intensities
                            ec_m_intensity=sum(ec_m_intensities)
                            ec_l_intensity=sum(ec_l_intensities)
                            ec_k_intensity=sum(ec_k_intensities)
                            bplus_intensity=sum(numpy.array(decay['beta+_intensity_list']))
                            tot_intensity=ec_intensity+bplus_intensity
                            if tot_intensity > 0.:
                                ec_intensity=ec_intensity/tot_intensity
                                ec_m_intensity=ec_m_intensity/tot_intensity
                                ec_l_intensity=ec_l_intensity/tot_intensity
                                ec_k_intensity=ec_k_intensity/tot_intensity
                                bplus_intensity=bplus_intensity/tot_intensity
                                ec_m_intensities=ec_m_intensities/tot_intensity
                                ec_l_intensities=ec_l_intensities/tot_intensity
                                ec_k_intensities=ec_k_intensities/tot_intensity
                                
                            else: #we presume a bplus decay
                                bplus_intensity=1.
                            decay['ec_cm+_intensity_list']=ec_m_intensities
                            decay['ec_cl_intensity_list']=ec_l_intensities
                            decay['ec_ck_intensity_list']=ec_k_intensities   
                            
                            possible_channels=[]
                            br_channels=[]
                            if bplus_intensity >0 :
                                possible_channels+=["BetaPlus"]
                                br_channels+=[br*bplus_intensity]
                            if ec_intensity >0 :
                                if (ec_m_intensity>0):
                                    possible_channels=["MshellEC"]+ possible_channels
                                    br_channels=[br*ec_m_intensity]+br_channels
                                if (ec_k_intensity>0):
                                    possible_channels+=["KshellEC"]
                                    br_channels+=[br*ec_k_intensity]
                                if (ec_l_intensity>0):
                                    possible_channels+=["LshellEC"]
                                    br_channels+=[br*ec_l_intensity]
                            
                            name_occuring_decay_channels+=possible_channels
                            br_channels_vec+=br_channels
                            for ll in range(len(possible_channels)):
                                occuring_decays_for_channels+=[decay]
                
            tot_br=numpy.sum(numpy.array(br_channels_vec))       
            if tot_br<.9999999:
                for name in decay_names:
                    name_wallet=name
                    if name=="BetaPlusEC":
                        name_wallet="BetaPlus"
                    if ( name not in name_occuring_decays and decay_wallet is not None and name_wallet in decay_wallet['type_list']):
                        k=np.where(np.array(decay_wallet['type_list'])==name_wallet) 
                        br=np.array(decay_wallet['BR_list'])[k]/100.
                        decay_ensdf_equivalent={"br":br,
                                            "Qval_ensdf":np.array(decay_wallet['Q_list'])[k]}
                        map_to_list_name={"BetaMinus":"betam_intensity_list", "BetaPlus":"beta+_intensity_list",
                                      "Alpha":"alpha_intensity_list"}
                    
                        decay_ensdf_equivalent[map_to_list_name[name_wallet]]=[1.]
                        decay_ensdf_equivalent["daughter_level_energy_list"]=[0.]
                        decay_ensdf_equivalent["g4beta_type_flag_list"]=[""]
                        name_occuring_decay_channels+=[name_wallet]
                        br_channels_vec+=[br]
                        occuring_decays_for_channels+=[decay_ensdf_equivalent]
                
               
                    
            #write the branching ratio of the different decay channel
            br_channels_vec=numpy.array(br_channels_vec)
            tot_br=numpy.sum(br_channels_vec)
            """
            print tot_br
            if tot_br !=1.:
                br_channels_vec= br_channels_vec/tot_br
            """
            k=0
            print name_occuring_decay_channels
        
            for channel in name_occuring_decay_channels:
                if channel=="ISO":
                    channel="IT"
                lines+=string_with_40_blanks[0:34-len(channel)]+channel+"       0.0000"+"   %.4e"%(br_channels_vec[k]/tot_br)+"\n"
                k+=1
            
            #write the decay intensity for the different channels    
            ll=0
            print "here ",name_occuring_decay_channels 
            for channel in name_occuring_decay_channels:
                decay=occuring_decays_for_channels[ll]
                print decay
                if (channel not in ["IT","ISO"]): #fort IT no need to write the intensities
                    Qval_ensdf=0.
                    if ("Qval_ensdf" in decay):
                        Qval_ensdf=decay["Qval_ensdf"]
                
                    Qval_qcalc=0.
                    if ("Qval_qcalc" in decay):
                        Qval_qcalc=decay["Qval_qcalc"]
                    elif ("Qval_EC_qcalc" in decay):
                        Qval_qcalc=decay["Qval_EC_qcalc"]
                        #For beta plus the Q-2emass is consider on G4 therefore  Qval EC should be taken
                        
                    
                    
                    daughter_level_energy_list=decay["daughter_level_energy_list"]
                    map_to_list_name={"BetaMinus":"betam_intensity_list", "BetaPlus":"beta+_intensity_list",
                                      "MshellEC":'ec_cm+_intensity_list',
                                      "LshellEC":'ec_cl_intensity_list',
                                      "KshellEC":'ec_ck_intensity_list',"Alpha":"alpha_intensity_list"}
                    intensity_list=numpy.array(decay[map_to_list_name[channel]])+0.000000000000000000000000001
                    intensity_list=100.*intensity_list*br_channels_vec[ll]/sum(intensity_list)
                    g4beta_type_flag_list=[]
                    if "Beta" in channel:
                        g4beta_type_flag_list=decay["g4beta_type_flag_list"]
                    print daughter_level_energy_list
                    for kk in range(len(daughter_level_energy_list)):
                        level_energy=daughter_level_energy_list[kk]
                        intensity=intensity_list[kk]
                        elevel_str="%.4f" %(level_energy)
                        elevel_str=string_with_20_blanks[0:13-len(elevel_str)]+elevel_str
                        intensity_str="%.4e" %(intensity)
                        intensity_str=string_with_20_blanks[0:13-len(intensity_str)]+intensity_str
                        g4beta_type_flag=""
                        if "Beta" in channel:
                            g4beta_type_flag=g4beta_type_flag_list[kk]
                        
                        qval_str="%.4f" %(Qval_qcalc+parent_elevel-level_energy)
                        qval_str="%.4f" %(Qval_ensdf+parent_elevel-level_energy)
                        """
                        if channel == "BetaPlus":
                            qval_str="%.4f" %(decay["Qval_beta+_qcalc"]+parent_elevel-level_energy)
                        """    
                        qval_str=string_with_20_blanks[0:13-len(qval_str)]+qval_str
                        if intensity>0 :
                            if ("Beta" in channel and g4beta_type_flag not in ["","allowed"]):
                                lines+=string_with_40_blanks[0:34-len(channel)]+channel+elevel_str+intensity_str+qval_str+"    "+g4beta_type_flag+"\n"
                            else:
                                lines+=string_with_40_blanks[0:34-len(channel)]+channel+elevel_str+intensity_str+qval_str+"\n"
                           
                    
                
                
                ll+=1
   
                
        
        
    if len(lines)>0:
        if not os.path.exists(new_directory):
            os.makedirs(new_directory)
        file_object=open(new_directory+"/z%i.a%i" %(Z,A),'w')
        file_object.write(lines)
        file_object.close()

## Write the file of the G4PhotoEvaporation database for a given nucleus
#
#  @param Z atomic/charge number, A mass number,
#  @param new_directory directpry where the new file with name z*.a* will be set,
#
def WriteG4PhotoEvaporationFile(Z,A,new_file_name,
                                      verbose=False,add_info_from_decay_file=True,parallel_run_nb=None):
    dic=ENSDF.ReadAdoptedLevelAndGammasInENSDFZipFile(Z,A,add_info_from_decay_file=add_info_from_decay_file)
    
    level_energies=numpy.array(dic["level_energies"])
    level_half_lifes=dic["half_life_list"]
    gammas=dic["gammas"]
    file_object=None

    for index_level in numpy.argsort(level_energies):
        level_energy=level_energies[index_level]
        gamma_energies=numpy.array(dic["gamma_energies"][index_level])
        gammas=dic["gammas"][index_level]
        level_half_life=level_half_lifes[index_level]
        JPi_vec_mother=dic["JPi_vec_list"][index_level]
        for index_gamma in numpy.argsort(gamma_energies):
            gamma_energy=gamma_energies[index_gamma]
            gamma=gammas[index_gamma]
            gamma_intensity=0.
            if ("RI" in gamma):
                gamma_intensity=gamma["RI"]
            daughter_levelE=level_energy-gamma_energy
            daughter_level_index=numpy.argmin(numpy.abs(level_energies-daughter_levelE))
            JPi_vec_daughter=dic["JPi_vec_list"][daughter_level_index]
            DJpi_str, MultiVec= GetDJpiAndPossibleMultipolarity(JPi_vec_mother,JPi_vec_daughter)
            if DJpi_str ==None:
                DJpi_str="99"
            J=0
            if (len(JPi_vec_mother)>0):
                J=JPi_vec_mother[0][0]
            M_gamma="-1"
            if "MP" in gamma:
                M_gamma=gamma["MP"]
            M_ratio=None
            if "MR" in gamma:
                M_ratio=gamma["MR"]
            #If the multipolarity is not given it takes the first of the computed MVec and put it in the datacard
            #before running BrIcc
            if M_gamma=="-1":
                M_ratio="-1"
                if (MultiVec is not None):
                    M_gamma="%s" %(MultiVec[0])
                    data_card=gamma["card"]
                    n0=32
                    data_card=data_card[0:n0]+M_gamma+data_card[n0+len(M_gamma):]
                    gamma["card"]=data_card
            str_jpi_in_file="%.2f" %(J+100)
            str_jpi_in_file=str_jpi_in_file[1:]
            if str_jpi_in_file[0]=="0":
                str_jpi_in_file=" "+str_jpi_in_file[1:]
            str_to_write="%.10e %.6e %.3e %s %.2e %s" %(level_energy,
                                                       gamma_energy,
                                                    gamma_intensity,
                                                        DJpi_str,
                                                        level_half_life,
                                                        str_jpi_in_file)
        
                
            total_conv_alpha=-1.
            if "CC" in gamma:
                total_conv_alpha=gamma["CC"]
                if (type(total_conv_alpha) == type({})):
                    total_conv_alpha=gamma["CC"]["value"]
            coeff_vec=numpy.zeros(10.)
            data_card=gamma["card"]
                    #print data_card
            icc_vec,ipf=Bricc.RunLocalBriccWithEnsdfFile(data_card,verbose=verbose,parallel_run_nb=parallel_run_nb)
            total_icc=numpy.sum(icc_vec)
            if (total_icc >0.):                                             
                coeff_vec=icc_vec/total_icc
            if total_conv_alpha <=0.:
                total_conv_alpha= total_icc                 
            str_to_write+=" %.3e" %(total_conv_alpha)
            for k in range(10):
                str_to_write+=" %.3e" %(coeff_vec[k])
            #str_to_write+=" %s %s" %(M_gamma,M_ratio)
            str_to_write+="\n"
            if file_object is None:
                file_object=open(new_file_name,'w')
            file_object.write(str_to_write)
    if file_object is not None:
        print new_file_name
        file_object.close() 

## Write the whole G4PhotoEvaporation database
#
#  @param [Zmin,Zmax] define the charge number range of nuclei for which files will be computed
#  @param new_directory directory where the new file with name z*.a* will be set,
#  @param parallel_run_nb integer number used to run BRICC in parallel
def WriteAllG4PhotoEvaporationDatabase(new_directory="../PhotoEvaporation3.0",
                                   Zmin=1,Zmax=100,parallel_run_nb=1):
   
    if not os.path.exists(new_directory):
            os.makedirs(new_directory)
    for Z in np.arange(Zmax-Zmin+1)+Zmin:
        for A in np.arange(3*Zmax-Zmin+1)+Zmin:
            if Z>=Zmin and Z<=Zmax:
                new_file_name=new_directory+"/z%i.a%i" %(Z,A)
                WriteG4PhotoEvaporationFile(Z,A,new_file_name,parallel_run_nb=parallel_run_nb)
   
## Write the whole G4RadioactiveDecay database
#
#   @param [Zmin,Zmax] define the charge number range of nuclei for which files will be computed
#   @param new_directory directory where the new file with name z*.a* will be set,
#           
def WriteAllG4RadDecayDatabase(new_directory="../RadioactiveDecay4.0",Zmin=1,Zmax=100):
    map_AZ_to_names, map_names_to_AZ,A_vec,Z_vec = NuclearWalletCards.ReadNucleusNameList()
    #list_decay_names,list_level_names=
    for i in range(len(Z_vec)):
        Z=Z_vec[i]
        A=A_vec[i]
        if Z>=Zmin and Z<=Zmax:
            WriteG4RadDecayFile(Z,A,new_directory=new_directory)

## Write all G4RadioactiveDecay files that have beta plus decay
#
#   @param [Zmin,Zmax] define the charge number range of nuclei for which files will be computed
#   @param new_directory directory where the new file with name z*.a* will be set,
#
def WriteAllG4RadDecayFilesWithBetaPlus(new_directory="../RadioactiveDecay4.0",Zmin=1,Zmax=20):
    map_AZ_to_names, map_names_to_AZ,A_vec,Z_vec = NuclearWalletCards.ReadNucleusNameList()
    for i in range(len(Z_vec)):
        Z=Z_vec[i]
        A=A_vec[i]
        if Z>=Zmin and Z<=Zmax:
            decays=ReadECBplusDecayInENSDFZipFile(Z,A)
            if decays is not None:
                print Z,A
                WriteG4RadDecayFile(Z,A,new_directory=new_directory)
   

#
#Functions to read and compute decays from the database
#

## Read radioactive decay data in a given block of G4RadioactiveDecay file
#
# This functions is used by ReadG4RadioactiveDecayFile
def ReadG4RadioactiveDecayBlock(lines):
    words=lines[0].split()
    decays={}
    decays['e_excitation']=float(words[1])
    decays['T1/2']=float(words[2])
    tot_intensity=0.
    tot_br_ratios=0.
    if len(lines) >1:
        for line in lines[1:]:
            if len(line) >13 and line[0]==" ":
                sub_line=line[13:]
                words=sub_line.split()
                type=words[0]
                if type not in decays:
                    decays[type]={}
                    decays[type]["branching_ratios"]=[]
                    decays[type]["daughter_ex"]=[]
                    decays[type]["Q"]=[]
                if (len(words)==3):
                    decays[type]["intensity"]=float(words[2])
                    tot_intensity+=float(words[2])
                if  (len(words)>3):
                    decays[type]["branching_ratios"]+=[float(words[2])]
                    decays[type]["daughter_ex"]+=[float(words[1])]
                    decays[type]["Q"]+=[float(words[3])]
                    tot_br_ratios+=float(words[2])
    decays["tot_intensity"]=tot_intensity
    decays["tot_br_ratios"]=tot_br_ratios
    return decays

## Read the radioactive decay data for a given nucleus into the G4RadioactiveDecay Database 
#
#   @param Z charge number of the nucleus
#   @param A mass number of the nucleus 
#   @param Z database_directory full path of the directory where the database is located 
def ReadG4RadioactiveDecayFile(Z,A,database_directory="../RadioactiveDecay4.0"):
    
    path="%s/z%i.a%i" %(database_directory,Z,A)
    if not os.path.exists(path):
        return None
    file_obj=open("%s/z%i.a%i" %(database_directory,Z,A),"r")
    lines=file_obj.readlines()
    file_obj.close()
    data_found=True
    n=len(lines)
    i=0
    blocks=None
    for line in lines:
        if line[0]=="P":
            if blocks is None:
                blocks=[[]]
            else:
                blocks+=[[]]
        if blocks is not None:
            blocks[-1]+=[line]
    
    decay_vec=[]
    
    if blocks is not None:
        for block in blocks:
            decay_vec+=[ReadG4RadioactiveDecayBlock(block)]
    else :
        return blocks
    
    return decay_vec

## Read the photo-evaporation data for a given nucleus into the G4PhotoEvaporation Database 
#
#   @param Z charge number of the nucleus
#   @param A mass number of the nucleus
#   @param database_directory  full path of the directory where the database is located 
def ReadG4PhotoEvaporationFile(Z,A,database_directory="../PhotoEvaporation3.0"):
    
    path="%s/z%i.a%i" %(database_directory,Z,A)
    if not os.path.exists(path):
        return None,None,None,None
    
    res=numpy.loadtxt(path,usecols=(0,1,2,6))
    e_levels=res[:,0]
    e_gammas=res[:,1]
    intensity_gammas=res[:,2]
    ic_alpha_factors=res[:,3]
    return e_levels,e_gammas,intensity_gammas,ic_alpha_factors

        

  
## Compute the decay products, with gamma deexitation, of a given nucleus contained in the G4 databases 
#
#  @param Z,A  atomic and mass number of  the nucleus
#  @param elevel energy level of the nucleus in keV (ground state as default)
#  @param database_directory  full path of the directory where the G4RadDecay database is located
#  @param photo_evap_database_directory full path of the directory where the photo-evaporation database is located
#           
def ComputeDecayProductsFromG4Database(Z,A,elevel=0.,database_directory="../g4data/RadioactiveDecay4.0",
                                        photo_evap_database_directory="../g4data/PhotonEvaporation3.0"):
    daughters={}  
    decay_vec=ReadG4RadioactiveDecayFile(Z,A,
                            database_directory=database_directory)
    alpha_energies=[]
    alpha_intensities=[]
    if decay_vec == None:
        decay_vec=[]
    for decay in decay_vec:
        if elevel == decay["e_excitation"]:
            print decay
            print decay.keys()
            for decay_type in ['MshellEC',"KshellEC","BetaPlus","LshellEC",
                               "BetaMinus","Alpha","IT"]:
                if decay_type in decay:
                    dZ,dA=GetDeltaZandAForDecayType(decay_type)
                    newZ=Z+dZ
                    newA=A+dA
                    print newZ,newA
                    daughter_id="z%i.a%i" %(newZ,newA)
                    if (daughter_id not in daughters):
                        daughters[daughter_id]={"e_levels":[],"intensities":[]}
                    daughter=daughters[daughter_id]
                    decay_channel=decay[decay_type]
                    print decay_channel
                    daughter_ex_vec=decay_channel['daughter_ex']
                        
                    br_vec=decay_channel['branching_ratios']
                    Q_vec=decay_channel["Q"]
                    for i in range(len(daughter_ex_vec)):
                        elevel_daughter=daughter_ex_vec[i]
                        br=br_vec[i]
                        Q=Q_vec[i]
                        
                        if elevel_daughter in daughter["e_levels"]:
                            ind=daughter["e_levels"].index(elevel_daughter)
                            daughter["intensities"][ind]+=br
                        else:
                            daughter["e_levels"]+=[elevel_daughter]
                            daughter["intensities"]+=[br]
                    
                        if decay_type=="Alpha":
                            mass_alpha=QCalc.GetMassFromQCALC(2,4,dir="../qcalc_data/august2012/",nuclear=True)/unit.keV
                            mass_parent=QCalc.GetMassFromQCALC(Z,A,dir="../qcalc_data/august2012/",nuclear=True)/unit.keV
                            mass_daughter=QCalc.GetMassFromQCALC(newZ,newA,dir="../qcalc_data/august2012/",nuclear=True)/unit.keV+elevel_daughter
                            #Q=mass_parent-mass_daughter-mass_alpha
                  
                            energy_alpha=(mass_daughter*Q+Q*Q/2.)/(Q+mass_daughter+mass_alpha)
                            print "alpha energy",energy_alpha
                            alpha_energies+=[energy_alpha]
                            alpha_intensities+=[br]
                            
    #No compute the gamma decays
    ###########################
    gamma_energies=[]
    gamma_intensities=[]
    for nuc_id in daughters:
        print nuc_id
        words=nuc_id.split(".")
        z=int(words[0][1:])
        a=int(words[1][1:])
        decay_vec1=ReadG4RadioactiveDecayFile(z,a,
                            database_directory=database_directory)
        e_levels_data,e_gammas_data,intensity_gammas_data,ic_alpha_factors_data=ReadG4PhotoEvaporationFile(z,a,database_directory=photo_evap_database_directory)
        if intensity_gammas_data is not None:
            intensity_gammas_data+=1.e-100
            
        daughter=daughters[nuc_id]
        if e_levels_data is not None:
            for i in range(len(daughter['e_levels'])):
                e_level=daughter['e_levels'][i]
                intensity=daughter['intensities'][i]
                e_gammas,intensity_gammas=ComputeGammaEmissionIntensitiesFromExcitedNucleus(
                        e_level,e_levels_data,e_gammas_data,intensity_gammas_data,ic_alpha_factors_data,decay_vec1)
            
                intensity_gammas=intensity_gammas*intensity
                for i in range(len(e_gammas)):
                    e=e_gammas[i]
                    intensity_gamma=intensity_gammas[i]
                    if e in gamma_energies:
                        index=gamma_energies.index(e)
                        gamma_intensities[index]+=intensity_gamma
                    else:
                        gamma_energies+=[e]
                        gamma_intensities+=[intensity_gamma]
    gamma_energies=np.array(gamma_energies)
    gamma_intensities=np.array(gamma_intensities)
    indices=np.argsort(gamma_energies)
    gamma_energies=gamma_energies[indices]
    gamma_intensities=gamma_intensities[indices]
    return daughters,gamma_energies,gamma_intensities,alpha_energies,alpha_intensities

## Function used internally to compute the gamma emission from an excited nucleus
#
#
def ComputeGammaEmissionIntensitiesFromExcitedNucleus(e_level,e_levels_data,e_gammas_data,
                                                      intensity_gammas_data,ic_alpha_factors_data,decay_vec):

    elevels_still_to_decay=[e_level]
    elevels_still_to_decay_intensity=[1.]
    intensity_gammas=[]
    e_gammas=[]
    
    elevels_in_decay_vec=[]
    if decay_vec is not None:
        for decay in decay_vec:
            elevels_in_decay_vec+=[decay['e_excitation']]
   
    intensity_electrons_ic_data=intensity_gammas_data*ic_alpha_factors_data
    print intensity_gammas
    while len(elevels_still_to_decay)>0:
        new_elevels_still_to_decay=[]
        new_elevels_still_to_decay_intensity=[]
        for k in range(len(elevels_still_to_decay)):
            e=elevels_still_to_decay[k]
            intensity=elevels_still_to_decay_intensity[k]
            k=0
            for elevel_in_decay in  elevels_in_decay_vec:
                print elevel_in_decay
                if np.abs(e-elevel_in_decay)<0.05:
                    factor=0.
                    if "IT" in decay_vec[k]:
                        factor=decay_vec[k]["IT"]["intensity"]
                    intensity*=factor
                k+=1
                    
                
            vec_for_search=np.abs(e_levels_data-e)
            arg_min=numpy.argmin(vec_for_search,axis=0)
            if (vec_for_search[arg_min]<0.5):
                indices=np.where(vec_for_search == vec_for_search[arg_min])
                
                selected_e=e_levels_data[indices]
                selected_gammas=e_gammas_data[indices]
                selected_intensity_gammas=intensity_gammas_data[indices]
                selected_alpha_factors=ic_alpha_factors_data[indices]
                selected_intensity_eic=intensity_electrons_ic_data[indices]
                relative_intensities=selected_intensity_gammas+selected_intensity_eic
                sum_intensities=np.sum(relative_intensities)
                selected_intensity_gammas=selected_intensity_gammas/sum_intensities
                relative_intensities=relative_intensities/sum_intensities
                for i in range(len(selected_gammas)):
                    e_gamma=selected_gammas[i]
                    intensity_gamma=intensity*selected_intensity_gammas[i]
                    e_daughter=e-e_gamma
                    intensity_daughter=relative_intensities[i]*intensity
                    if e_daughter>0.001:
                        index=-1
                        if e_daughter in new_elevels_still_to_decay:
                             index=new_elevels_still_to_decay.index(e_daughter)
                             new_elevels_still_to_decay_intensity[index]+=intensity_daughter
                        else:
                             new_elevels_still_to_decay+=[e_daughter]
                             new_elevels_still_to_decay_intensity+=[intensity_daughter]
                    if e_gamma in e_gammas:
                        index=e_gammas.index(e_gamma)
                        intensity_gammas[index]+=intensity_gamma
                    else:
                        e_gammas+=[e_gamma]
                        intensity_gammas+=[intensity_gamma]
        elevels_still_to_decay=new_elevels_still_to_decay
        elevels_still_to_decay_intensity=new_elevels_still_to_decay_intensity
                
    return np.array(e_gammas),np.array(intensity_gammas)
    
## Compute the full radioactive chain of a given nucleus
#
#  @param Z,A atomic and mass number of  the nucleus
#  @param elevel energy level of the nucleus in keV (ground state as default)
#  @param database_directory  full path of the directory where the G4RadDecay database is located
#  @param  photo_evap_database_directory full path of the directory where the photo-evaporation database is located 
def ComputeFullChain(Z,A,database_directory="../g4data/RadioactiveDecay4.0",
                                        photo_evap_database_directory="../g4data/PhotonEvaporation3.0"):
    daughters,gamma_energies,gamma_intensities,alpha_energies,alpha_intensities=ComputeDecayProductsFromG4Database(Z,A,0.,database_directory,
                                        photo_evap_database_directory)
    list_daughters=["z%i.a%i" %(Z,A)]
    while len(daughters)>0:
        list_new_daughters=[]
        for nuc_id in daughters:
            if nuc_id not in list_daughters:
                list_daughters+=[nuc_id]
            Z1=int(nuc_id.split(".")[0][1:])
            A1=int(nuc_id.split(".")[1][1:])
            daughters1,gamma_energies,gamma_intensities,alpha_energies,alpha_intensities=ComputeDecayProductsFromG4Database(Z1,A1,0.,database_directory,
                                        photo_evap_database_directory)
            for nuc_id1 in  daughters1:
                if nuc_id1 not in list_new_daughters and nuc_id1 not in daughters and nuc_id1 not in list_daughters:
                    list_new_daughters+=[nuc_id1]
                    list_daughters+=[nuc_id1]
        daughters=list_new_daughters
    for daughter in list_daughters[::-1]:
        ViewFile("%s/%s" %(database_directory,daughter))
        




                

def GetBetaInfoFromLines(lines):
    parent_elevel=0.
    bminus_parent_elevels=[]
    bminus_elevels=[]
    bminus_index_lines=[]
    
    bplus_parent_elevels=[]
    bplus_elevels=[]
    bplus_index_lines=[]
    
    for i in range(len(lines)):
        words=lines[i].split()
        if words[0]=="P":
            parent_elevel=float(words[1])
        if len(words)>3:
            if words[0] == "BetaMinus":
                bminus_parent_elevels+=[parent_elevel]
                bminus_elevels+=[float(words[1])]
                bminus_index_lines+=[i]
            if words[0] == "BetaPlus":
                bplus_parent_elevels+=[parent_elevel]
                bplus_elevels+=[float(words[1])]
                bplus_index_lines+=[i]
    return bminus_parent_elevels,bminus_elevels,bminus_index_lines,bplus_parent_elevels,bplus_elevels,bplus_index_lines
                

def AddForbiddenStatesToOldDatabase(new_database_dir="../RadioactiveDecay4.0",
                                    old_database_dir="/Users/laurent/g4data/RadioactiveDecay3.5",
                                    replace_data_base_dir="../RadioactiveDecay3.5_withForbiddenInfo"):
    #Read the forbidden states file
    obj_file=open(new_database_dir+"/Forbidden.txt",'r')
    lines=obj_file.readlines()
    old_file_obj=None
    old_lines=None
    old_file_name=""
    old_beta_plus_parent_elevel=[]
    old_beta_plus_elevels=[]
    old_beta_minus_elevels=[]
    
    old_beta_minus_parent_elevel=[]
    old_beta_plus_iline=[]
    old_beta_minus_iline=[]
    parent_level=-10000.
    old_parent_level=-5000000.
    list_cases_to_update=[]
    list_files_to_update=[]
    
    for line in lines:
        
        words=line.split()
        file_name=words[0][0:-1]
        intensity=float(words[3])
        if file_name not in list_files_to_update:
            list_files_to_update+=[file_name]
    
    for file_name in list_files_to_update:
        old_path=old_database_dir+"/"+file_name
        new_path=new_database_dir+"/"+file_name
        replace_path=replace_data_base_dir+"/"+file_name
        if os.path.exists(old_path):
            old_obj_file=open(old_path,'r')
            old_lines=old_obj_file.readlines()
            old_obj_file.close()
            
            bm_p_elevs,bm_elevs,bm_ilines,bp_p_elevs,bp_elevs,bp_ilines=GetBetaInfoFromLines(old_lines)
            bm_p_elevs=np.array(bm_p_elevs)
            bm_elevs=np.array(bm_elevs)
            bm_ilines=np.array(bm_ilines)
            bp_p_elevs=np.array(bp_p_elevs)
            bp_elevs=np.array(bp_elevs)
            bp_ilines=np.array(bp_ilines)
            
            
            new_obj_file=open(new_path,'r')
            new_lines=new_obj_file.readlines()
            new_obj_file.close()
            parent_elevel=0.
            for line in new_lines:
                words=line.split()
                if words[0]=="P":
                    parent_elevel=float(words[1])
                if words[-1] in ["firstForbidden", "uniqueFirstForbidden",
                                 "uniqueSecondForbidden","uniqueThirdForbidden"]:
                    if words[0] == "BetaMinus":
                        elev=float(words[1])
                        minimum_p=np.min(np.abs(bm_p_elevs-parent_elevel))
                        if minimum_p<1.:
                            indices=np.where(np.abs(bm_p_elevs-parent_elevel)==minimum_p)
                            elevs=bm_elevs[indices]
                            ilines=bm_ilines[indices]
                            minimum_e=np.min(np.abs(elevs-elev))
                            if minimum_e<.5:
                                ind=np.argmin(np.abs(elevs-elev))
                                i=ilines[ind]
                                new_line=old_lines[i][0:73]+line[73:]
                                old_lines[i]=new_line
                    if words[0] == "BetaPlus" and len(bp_p_elevs)>0:
                        elev=float(words[1])
                        minimum_p=np.min(np.abs(bp_p_elevs-parent_elevel))
                        if minimum_p<1.:
                            indices=np.where(np.abs(bp_p_elevs-parent_elevel)==minimum_p)
                            elevs=bp_elevs[indices]
                            ilines=bp_ilines[indices]
                            minimum_e=np.min(np.abs(elevs-elev))
                            if minimum_e<.5:
                                ind=np.argmin(np.abs(elevs-elev))
                                i=ilines[ind]
                                new_line=old_lines[i][0:73]+line[73:]
                                old_lines[i]=new_line
        
            file_obj=open(replace_path,'w')
            for old_line in old_lines:
                file_obj.write(old_line)
            file_obj.close()
                                
