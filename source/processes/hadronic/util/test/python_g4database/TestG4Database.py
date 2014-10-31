## @package TestG4Database
#  TestG4Database PYTHON module to test the  G4 RadDEcay and PhotoEvaporation database with other data
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

from QCalc import *
from NuclearWalletCards import *
from utilities import *
from Bricc import *
from ENSDF import *
from G4Database import *

uma_in_MeV=931.494028

## Compare G4RadDecay database and nuclear wallet data for a given nucleus 
#
#
def CompareG4AndNuclearWalletFile(Z,A,nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../RadioactiveDecay4.0",
                                 photo_evap_database_directory="../g4data/PhotonEvaporation3.0"):
    list_stable_nuclei, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic=ReadNuclearWalletFile(name_file=nuc_wallet_file)
    G4decay_vec=ReadG4RadioactiveDecayFile(Z,A,database_directory="../RadioactiveDecay4.0")
    name="z%i.a%i" %(Z,A)
    is_in_G4_raddecay=len(G4decay_vec) >0
    is_in_NuclearWallet=name in list_radioactive_nuclei
    are_the_same=False
    res_comp=[]
    identical=True
    if is_in_G4_raddecay and is_in_NuclearWallet:
        wallet_decays=radioactive_nuc_dic[name]
        print wallet_decays
        excited_energy_vec=np.array(wallet_decays["excited_energy_vec"])
        for G4decay in G4decay_vec:
            print G4decay
            res_comp+={}
            energy=G4decay['e_excitation']
            indices=np.where(np.abs(excited_energy_vec-energy)<=0.001)
            if len(indices[0]) ==1:
                e=excited_energy_vec[indices][0]
                WalletDecay=wallet_decays[e]
                print WalletDecay
            else:
                identical=False
## Compare G4RadDecay database and nudat2 data for a given nucleus 
#
#
def CompareG4DatabaseOutputAndNudat2Output(Z,A,nuc_name,elevel=0.,g4rad_decay_dir="../g4data/RadioactiveDecay4.0",
                                           photo_evap_database_directory="../g4data/PhotonEvaporation3.0",show=False):
    
    daughters,gamma_energies,gamma_intensities,alpha_energies,alpha_intensities=ComputeDecayProductsFromG4Database(Z,A,elevel=elevel,
                                                                                    database_directory=g4rad_decay_dir,
                                                                              photo_evap_database_directory=photo_evap_database_directory)

    #print np.sum(daughters['z17.a34']['intensities'])
    pl.figure(figsize=(6,6))
    pl.rcParams['legend.fontsize']=10
    pl.subplot(2,1,1)

    indices=np.where(gamma_intensities>1.e-3)
    gamma_energies=gamma_energies[indices]
    gamma_intensities=gamma_intensities[indices]
    
    min_gamma=0
    max_gamma=1
    
    if len(gamma_energies)>0:
        min_gamma=np.min(gamma_energies)*0.9
        max_gamma=np.max(gamma_energies)*1.1
        
    
    Energy_vec_g,Intensity_vec_g=GetGammaIntensitiesFromRadDecayNudat2(Z,A,dir="../nudat2_data/nuclear_decay/november2013",with_xrays=False)
    
    
    for i in range(len(gamma_energies)):
        print "e_gamma  G4 ",gamma_energies[i], "intensity ",gamma_intensities[i]
        e=gamma_energies[i]
        id=np.argmin(np.abs(Energy_vec_g-e))
        print "Nudat2 e_gamma ",Energy_vec_g[id]," intensity ",Intensity_vec_g[id]
        intensity=gamma_intensities[i] 
        if i==0:
            pl.semilogy([e,e],[1.e-5,intensity],"k",label="G4Radecay4.0\nG4PhotonEvaporation3.0")
        else:
            pl.semilogy([e,e],[1.e-5,intensity],"k")
            

    pl.subplot(2,1,2)
    Energy_vec,Intensity_vec=GetAlphaIntensitiesFromRadDecayNudat2(Z,A,dir="../nudat2_data/nuclear_decay/august2012")
    print Energy_vec,Intensity_vec
    for i in range(len(alpha_energies)):
        e=alpha_energies[i]
        intensity=alpha_intensities[i] 
        print "G4 e_alpha ",alpha_energies[i]," intensity ",alpha_intensities[i] 
        if len(Energy_vec):
            id=np.argmin(np.abs(Energy_vec-e))
            print "Nudat2 e_alpha ",Energy_vec[id]," intensity ",Intensity_vec[id]
        pl.semilogy([e,e],[1.e-5,intensity],"k")
    
    

    
    
    
    
    pl.semilogy(Energy_vec_g,Intensity_vec_g,"b*",label="nudat2")
    pl.ylim([1.e-3,100.])
    
    print min_gamma,max_gamma
    pl.subplot(2,1,1)
    pl.xlim([min_gamma,max_gamma])
    pl.xlabel("Ekin [keV]")
    pl.ylabel("Intensity [%]")
    pl.title("Gamma intensity from %s decay" %(nuc_name))
    pl.legend()

    pl.subplot(2,1,2)
   
    Energy_vec,Intensity_vec=GetAlphaIntensitiesFromRadDecayNudat2(Z,A,dir="../nudat2_data/nuclear_decay/august2012")
    
    pl.semilogy(Energy_vec,Intensity_vec,"b*")
    
    pl.ylim([1.e-3,100.])
    if show:
        pl.show()
    if not os.path.exists("../plots/ValidationG4RadDatabase4.1_PhotoEvap3.0"):
        os.mkdir("../plots/ValidationG4RadDatabase4.1_PhotoEvap3.0")
    
    pl.savefig("../plots/ValidationG4RadDatabase4.1_PhotoEvap3.0/GammaAndAlphaProduction%s.pdf" %(nuc_name) ) 
    
    
## Find missing nuclei in G4RadDEcay database compare to nuclear wallet card
#
#                
def FindMissingNucleiInG4Rdecay(nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../RadioactiveDecay4.0"):
    list_stable_nuclei, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic=ReadNuclearWalletFile(name_file=nuc_wallet_file)
    g4decay_files=os.listdir(g4rad_decay_dir)
    list_missing_nuclei=[]
    list_additional_nuclei=[]
    for nucleus_name in list_radioactive_nuclei:
        if nucleus_name not in g4decay_files:
            dic=radioactive_nuc_dic[nucleus_name]
            is_g4rad_candidate=False
            for e in dic['excited_energy_vec']:
                for dec_mode in ['B-','EC','A','IT']:
                    if dec_mode in dic[e]:
                        if dic[e][dec_mode]>0.:
                            is_g4rad_candidate=True
            if is_g4rad_candidate:
                list_missing_nuclei+=[nucleus_name]
    for nucleus_name in g4decay_files:
        if nucleus_name not in list_radioactive_nuclei:
            list_additional_nuclei+=[nucleus_name]
    return list_missing_nuclei

## Compare List of files in two different versions of G4RadDecay database
#
#
def CompareListFilesInG4RaddecayDatabases(g4rad_decay_dir1="../RadioactiveDecay4.0",
                                          g4rad_decay_dir2="../g4data/RadioactiveDecay3.7"):
    g4decay_files1=os.listdir(g4rad_decay_dir1)
    g4decay_files2=os.listdir(g4rad_decay_dir2)
    list_common_file=[]
    list_extra_files_in_database1=[]
    list_extra_files_in_database2=[]
    for nucleus_name in g4decay_files1:
        if nucleus_name not in g4decay_files2:
            list_extra_files_in_database1+=[nucleus_name]
        else:
            list_common_file+=[nucleus_name]
    for nucleus_name in g4decay_files2:
        if nucleus_name not in g4decay_files1:
            list_extra_files_in_database2+=[nucleus_name]
    return list_extra_files_in_database1,list_extra_files_in_database2
   
            
## Compare two alpha decays from tow different versions of the database for a given nucleus
#
#
def CompareG4RadioactiveDecayFilesForAlpha(Z,A,database_directory="../RadioactiveDecay4.0",
                                           database_directory1="../final_g4data/RadioactiveDecay4.0"):
    decay_vec=ReadG4RadioactiveDecayFile(Z,A,database_directory)
    print "vec1 ",decay_vec
    decay_vec1=ReadG4RadioactiveDecayFile(Z,A,database_directory1)
    print "vec2 ",decay_vec1
    
## Check if  a file of the radioactive decay database gives empty branching ratio 
#
#
def CheckEmptyBranchingRatioG4RadioactiveDecayFilesForAlpha(Z,A,
                                           database_directory="../final_g4data/RadioactiveDecay4.0"):
    decay_vec=ReadG4RadioactiveDecayFile(Z,A,database_directory)
    tot_br_ratios_is_ok=True
    is_alpha=False
    with_decay=False
    if decay_vec is not None:
        with_decay=True
        for decay in decay_vec:
            if 'Alpha' in decay:
                is_alpha=True
                if len(decay['Alpha']['branching_ratios']) == 0:
                    tot_br_ratios_is_ok=False
    return with_decay,is_alpha,tot_br_ratios_is_ok


## Find files of the radioactive decay database with empty branching ratio for alpha decays
#
#
def CheckAllG4RadDecayFilesForAlpha(database_directory="../final_g4data/RadioactiveDecay4.0_patch1",Zmin=1,Zmax=100):
    map_AZ_to_names, map_names_to_AZ,A_vec,Z_vec = ReadNucleusNameListFromNucWalletNudat2()
    #list_decay_names,list_level_names=
    for i in range(len(Z_vec)):
        Z=Z_vec[i]
        A=A_vec[i]
        if Z>=Zmin and Z<=Zmax:
            
             with_decay,is_alpha,tot_br_ratios_is_ok=CheckEmptyBranchingRatioG4RadioactiveDecayFilesForAlpha(Z,A,
                                           database_directory)
             if not tot_br_ratios_is_ok:
                print Z,A
                """
                WriteG4RadDecayFile(Z,A,new_directory="../final_g4data/RadioactiveDecay4.1_patch1",
                        zip_file_name_decays="../ensdf_data/september_2013/decays.zip",
                                     zip_file_name_levels="../ensdf_data/august_2012/all_levels.zip")
                """
                
## Find files with missing info in the old databases
#
#           
def FindFilesWithMissingInfoInOldDatabase(new_database_dir="../RadioactiveDecay4.0",
                                    old_database_dir="/Users/laurent/g4data/RadioactiveDecay3.6",
                                    store_directory="../temp"):
    files=os.listdir(old_database_dir)
    list_files=[]
    for file in files:
        file_obj=open("%s/%s" %(old_database_dir,file),"r")
        lines_old=file_obj.readlines()
        file_obj.close()
        last_line=lines_old[-1]
        if last_line[0]=="P":
            #print file
            list_files+=[file]
            #print last_line
            file_store_obj=open("%s/%s_old" %(store_directory,file),"w")
            for line in lines_old:
                file_store_obj.write(line)
            file_store_obj.close()
            if os.path.exists("%s/%s" %(new_database_dir,file)):
                file_obj1=open("%s/%s" %(new_database_dir,file),"r")
                lines=file_obj1.readlines()
                file_store_obj=open("%s/%s_new" %(store_directory,file),"w")
                
                last_P_index=-1
                i=0
                for line in lines:
                    file_store_obj.write(line)
                    if line[0]=="P":
                        last_P_index=i
                    i+=1
                #print lines[last_P_index]
                file_store_obj.close()
            
            #print last_line[0:-1]
   
## Find files with some missing info 
#
#
def CheckPotentialFilesWithMissingInfo(database_dir="../RadioactiveDecay4.0",
                                    store_directory="../temp"):
    files=os.listdir(database_dir)
    list_files=[]
    for file in files:
        file_obj=open("%s/%s" %(database_dir,file),"r")
        lines=file_obj.readlines()
        file_obj.close()
        last_line=lines[-1]
        is_bad_file=False
        data_found=True
        n=len(lines)
        i=0
        while (i<n and not is_bad_file):
            line=lines[i]
            i+=1
            if line[0]=="P":
                if not data_found:
                    is_bad_file=True
                else:
                    data_found=False
            else:
                if line[0]==" " and len(line.split())>=3:
                    data_found=True
        if not data_found:
            is_bad_file=True
        if is_bad_file:
            list_files+=[file]
    return list_files
## Find files without  decay info 
#
#                
def CheckPotentialFilesWithoutDecayInfo(database_dir="../RadioactiveDecay4.0",
                                    store_directory="../temp"):
    files=os.listdir(database_dir)
    list_files=[]
    for file in files:
        file_obj=open("%s/%s" %(database_dir,file),"r")
        lines=file_obj.readlines()
        file_obj.close()
        last_line=lines[-1]
        is_bad_file=False
        data_found=False
        n=len(lines)
        i=0
        while (i<n and data_found==False):
            line=lines[i]
            i+=1
            if line[0]==" " and len(line.split())>=3:
                data_found=True
        if not data_found:
            list_files+=[file]
    return list_files                 
                    


if __name__ == "__main__":
    x=1

    #WriteG4RadDecayFile(99,250)
    """
    WriteG4RadDecayFile(38,87)
    WriteG4RadDecayFile(35,80)
    WriteG4RadDecayFile(19,40)
    """ 
#    WriteAllG4RadDecayFile(Zmin=21)
    
#    WriteG4RadDecayFile(21,56)
    """
    WriteG4RadDecayFile(89,228)
    decay_list=ReadECBplusDecayInENSDFZipFile(89,228)
    print     decay_list
    """  
    """
    dic=ReadAdoptedLevelAndGammasInENSDFZipFile(91,234,zip_file_name_levels="../ensdf_data/august_2012/all_levels.zip",
                                                     zip_file_name_decays="../ensdf_data/august_2012/all_decays.zip",
                                                     add_info_from_decay_file=True)
    print dic["gammas"]
    print dic.keys()

    for i in range(len(dic["ensdf_level_energies"])):
        ref_tag=dic["level_reference_list"][i]
        e_ref=dic["ReferenceLevels"][ref_tag]
        
        print "Level",dic["ensdf_level_energies"][i],ref_tag,dic["ensdf_level_energies"][i]+e_ref
        
        print dic["gamma_energies"][i]
        print dic["gammas"][i]
    print dic["ReferenceLevels"]
    """
    
    
    #WriteG4RadDecayFile(10,24,new_directory="../RadioactiveDecayTest")
    """
    decay_vec=ReadG4RadioactiveDecayFile(10,24)
    print decay_vec
    decay_vec=ReadG4RadioactiveDecayFile(9,24)
    print decay_vec
    list_stable_nuc, list_radioactive_nuclei, radioactive_nuc_dic, stable_nuc_dic=ReadNuclearWalletFile()
    print list_stable_nuc
    print list_radioactive_nuclei
    """
    """
    FindMissingNucleiInG4Rdecay(nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../g4data/RadioactiveDecay3.7")
    
    CompareListFilesInG4RaddecayDatabases(g4rad_decay_dir1="../RadioactiveDecay4.0",
                                          g4rad_decay_dir2="../g4data/RadioactiveDecay3.7")
    
    nuc_name, allowed_decays, decays_dic, Parent_elevel_str_list, Parent_elevel_list, \
                           Parent_half_life_str_list, Parent_half_life_list = ReadDecaysInENSDFZipFile(84,214)
    
    decay=ReadAlphaDecayInENSDFZipFile(84,214,zip_file_name="../ensdf_data/august_2012/all_decays.zip")
    print decay
    
    """

    #UpdateDatabaseWithMissingIsomers(new_directory="../RadioactiveDecay4.0")
    """
    file_names=["z11.a25","z27.a56","z53.a131","z54.a131","z81.a206","z82.a209",
                "z87.a221","z88.a225","z89.a227","z91.a231","z91.a233","z12.a25","z26.a57",
                "z9.a19","z9.a21","z11.a22","z12.a29","z13.a28","z13.a29","z12.a31"]
    """
    """
    file_names=["z18.a37"]
    for file_name in file_names:
        print file_name
        ViewFile("../RadioactiveDecay4.0/%s" %(file_name))
        print "**********************"
        ViewFile("../g4data/RadioactiveDecay3.7/%s" %(file_name))
        print "**********************"
    nuc_name, allowed_decays, decays_dic, Parent_elevel_str_list, Parent_elevel_list, \
                    Parent_half_life_str_list, Parent_half_life_list \
                                                    =ReadDecaysInENSDFZipFile(18,37)
    """
            
    """
    """
    """
    CompareG4AndNuclearWalletFile(9,24,nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../RadioactiveDecay4.0")
    """
    #WriteG4RadDecayFile(9,24,new_directory="../RadioactiveDecay4.0")
    #WriteG4RadDecayFile(84,214,new_directory="../RadioactiveDecay4.0")
    """
    AddItemInListMissingIsomers(83,210,46.539)
    AddItemInListMissingIsomers(92,234,1552.6200)
    print ids_isomer_to_add_in_database
    print elevels_isomer_to_add_in_database
    print elevels_daughter_isomer_to_add_in_database
    print half_life_isomer_to_add_in_database
    WriteG4RadDecayFile(83,210,new_directory="../RadioactiveDecay4.0")
    WriteG4RadDecayFile(92,234,new_directory="../RadioactiveDecay4.0")
    """
    #print allowed_decays
    
    #WriteG4RadDecayFile(23,43,new_directory="../RadioactiveDecay4.0")
    #WriteG4RadDecayFile(2,5,new_directory="../RadioactiveDecay4.0")
    #WriteG4RadDecayFile(74,186,new_directory="../RadioactiveDecay4.0")
    """
    print FindMissingNucleiInG4Rdecay(nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../RadioactiveDecay4.0")
    
    print FindMissingNucleiInG4Rdecay(nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../g4data/RadioactiveDecay3.7")
    """
    
    #WriteAllG4RadDecayFileWithBetaPlus(new_directory="../RadioactiveDecay4.0",Zmin=69,Zmax=100)
    
    """
    #ReadBminusDecayInENSDFZipFile(10,24,zip_file_name="../ensdf_data/august_2012/all_decays.zip")
    
    WriteG4RadDecayFile(90,234,new_directory="../NewRadioactiveDecayAugust2012")
    WriteG4RadDecayFile(91,234,new_directory="../NewRadioactiveDecayAugust2012")

    WriteG4RadDecayFile(92,235,new_directory="../NewRadioactiveDecayAugust2012")
    
    WriteG4PhotoEvaporationFile(90,231,"../new_photoevaporation_august2012/z90.a231",
                                      verbose=False,add_info_from_decay_file=True)
    
    WriteG4PhotoEvaporationFile(90,234,"../new_photoevaporation_august2012/z90.a234",
                                      verbose=False,add_info_from_decay_file=True)
    
    WriteG4PhotoEvaporationFile(64,158,"../new_photoevaporation_august2012/z64.a158",
                                      verbose=False,add_info_from_decay_file=True)
    
    WriteG4PhotoEvaporationFile(12,24,"../PhotoEvaporation3.0/z12.a24",
                                      verbose=False,add_info_from_decay_file=True)
    
    
    WriteG4PhotoEvaporationFile(48,111,"../PhotoEvaporation3.0/z48.a111",
                                      verbose=False,add_info_from_decay_file=True)
    
    WriteAllG4PhotoEvaporationFile(new_directory="../PhotoEvaporation3.0",
                                   Zmin=61,Zmax=100)
    
    
    WriteAllG4PhotoEvaporationFile(new_directory="../new_photoevaporation_august2012",
                                   Zmin=1,Zmax=100)
    
    WriteG4RadDecayFile(31,73,new_directory="../NewRadioactiveDecayAugust2012")
    WriteG4RadDecayFile(32,73,new_directory="../NewRadioactiveDecayAugust2012")
    WriteG4PhotoEvaporationFile(32,73,"../new_photoevaporation_august2012/z32.a73",
                                      verbose=False,add_info_from_decay_file=True)
    data_card_ensdf="231TH  G 543.66    3       32 7 M1                      0.178\n"
    
    icc_vec,a_vec=RunLocalBriccWithEnsdfFile(data_card_ensdf,verbose=False)
    print icc_vec 
    print GetValuesForJPiField("(3/2+,5/2-")   
    print GetValuesForJPiField("(3/2+ TO 7/2-)")

    AddForbiddenStatesToOldDatabase(new_database_dir="../RadioactiveDecay4.0",
                                    old_database_dir="/Users/laurent/g4data/RadioactiveDecay3.5")
    
    
    FindFilesWithMissingInfoInOldDatabase(new_database_dir="../RadioactiveDecay4.0",
                                    old_database_dir="/Users/laurent/g4data/RadioactiveDecay3.6")
    """
    
    #WriteG4RadDecayFile(2,7,new_directory="../RadioactiveDecay4.0")
    #print ReadDecaysInENSDFZipFile(100,242)
    """
    WriteG4PhotoEvaporationFile(16,80,"./z16.a36",
                                     verbose=False,add_info_from_decay_file=True)
    """
    #print GetNucleusName(16,80)
    #WriteG4RadDecayFile(16,80,new_directory="./")
    #WriteListOfITDecays()
    
    """
    WriteG4PhotoEvaporationFile(88,226,"../PhotoEvaporation3.0/z88.a230",
                                     verbose=False,add_info_from_decay_file=True)
    """
    Z_vec=[27,27,56,55,55,63,48,58,50,39,38,19,81,82,82,83,88,89,90,95,92,92,92,94,94,56,31,64,81,
       19,13,33,91,90,91,89,90,86,84,63]
    A_vec=[60,57,133,134,137,152,109,139,113,88,90,40,208,210,214,214,226,228,234,241,234,235,238,
       238,239,133,73,158,204,40,26,76,234,231,231,227,230,222,218,141]
    Nucname_vec=["Co60","Co57","Ba133","Cs134","Cs137",
             "Eu152","Cd109","Ce139","Sn113","Y88",
             "Sr90","K40","Tl208","Pb210","Pb214","Bi214","Ra226","Ac228","Th234","Am241","U234",
             "U235","U238","Pu238","Pu239","Ba133","Ga73","Gd158","204TL","K40","Al26","As76",
             "Pa234","Th231","Pa231","Ac227","Th230","Rn222","Po218","Eu141"
             ]

    ind_map={}
    i=0
    #Nucname_vec=["Eu152","Cd109","Ce139","Sn113","Y88"]
    for name in Nucname_vec:
        print name
        ind_map[name]=i
        i+=1
        
   
    """    
    for nuc_name in Nucname_vec:  
        i=ind_map[nuc_name]
        Z=Z_vec[i]
        A=A_vec[i]
        print i,Z,A
        ViewFile("../g4data/RadioactiveDecay3.7/z%i.a%i" %(Z,A))
        ViewFile("../g4data/RadioactiveDecay4.0/z%i.a%i" %(Z,A))
        print i,Z,A
        print nuc_name
        elevel=0.
        if nuc_name=="Pa234":
            elevel=73.29
        if Z<41 or Z>50:
            CompareG4DatabaseOutputAndNudat2Output(Z,A,nuc_name,elevel=elevel,
                                                   g4rad_decay_dir="../g4data/RadioactiveDecay4.0",
                                                   photo_evap_database_directory="../g4data/PhotonEvaporation3.0",show=False)
        ViewFile("../g4data/RadioactiveDecay4.0/z%i.a%i" %(Z,A))
        ViewFile("../g4data/RadioactiveDecay3.4/z%i.a%i" %(Z,A))
    
    Z=72
    A=163
    ViewFile("../g4data/PhotonEvaporation2.3/z%i.a%i" %(Z,A))
    ViewFile("../g4data/PhotonEvaporation3.0/z%i.a%i" %(Z,A))
    
    WriteG4PhotoEvaporationFile(72,163,"test",
                                     verbose=False,add_info_from_decay_file=True)
    ViewFile("test")
    """
       
    """
    Z=41
    A=92
    #WriteG4RadDecayFile(Z,A,new_directory="../g4data/RadioactiveDecay4.0")
    ViewFile("../g4data/RadioactiveDecay3.7/z%i.a%i" %(Z,A))
    ViewFile("../g4data/RadioactiveDecay4.0/z%i.a%i" %(Z,A))
    """
    """
    print ReadDecaysInWalletCards(Z,A)
    
    FindFilesWithMissingInfoInOldDatabase(new_database_dir="../g4data/PhotonEvaporation3.0",
                                    old_database_dir="../g4data/PhotonEvaporation2.3")
    
    list1,list2=CompareListFilesInG4RaddecayDatabases(g4rad_decay_dir1="../g4data/PhotonEvaporation3.0",
                                          g4rad_decay_dir2="../g4data/PhotonEvaporation2.3")
    for file in list2:
        print file
    
    list1,list2=CompareListFilesInG4RaddecayDatabases(g4rad_decay_dir1="../g4data/RadioactiveDecay4.0",
                                          g4rad_decay_dir2="../g4data/RadioactiveDecay3.7")
    for file in list2:
        print file
        ViewFile("../g4data/RadioactiveDecay4.0/%s" %(file))
        ViewFile("../g4data/RadioactiveDecay3.7/%s" %(file))
    
    list=FindMissingNucleiInG4Rdecay(nuc_wallet_file="../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt",
                                g4rad_decay_dir="../g4data/RadioactiveDecay4.0")
    for file in list:
        print file
    """
    Z=92
    A=222
    WriteG4RadDecayFile(Z,A,new_directory="../RadioactiveDecayTest",
                        zip_file_name_decays="../ensdf_data/september_2013/decays.zip",
                                     zip_file_name_levels="../ensdf_data/august_2012/all_levels.zip")
    CompareG4RadioactiveDecayFilesForAlpha(94,230,database_directory="../RadioactiveDecay4.0",
                                           database_directory1="../final_g4data/RadioactiveDecay4.0")
    print CheckEmptyBranchingRatioG4RadioactiveDecayFilesForAlpha(94,230,
                                           database_directory="../RadioactiveDecay4.0")
    
    Z=94
    A=230
    ViewFile(file_name="../RadioactiveDecay4.0/z%i.a%i" %(Z,A))
    ViewFile(file_name="../final_g4data/RadioactiveDecay4.1_patch1/z%i.a%i" %(Z,A))
    
    CheckAllG4RadDecayFilesForAlpha(database_directory="../final_g4data/RadioactiveDecay4.1_patch1",Zmin=1,Zmax=100)
    print "test"
    CheckAllG4RadDecayFilesForAlpha(database_directory="../final_g4data/RadioactiveDecay4.0_patch1",Zmin=1,Zmax=100)
    
    """
    CompareG4DatabaseOutputAndNudat2Output(96,237,"237Cm",elevel=0.,g4rad_decay_dir="../final_g4data/RadioactiveDecay4.1_patch1",
                                           photo_evap_database_directory="../g4data/PhotonEvaporation3.0",show=False)
    
    
    """
    """
    Z=90
    A=221
    WriteG4PhotoEvaporationFile(Z,A,"../PhotoEvaporationTest/Z%i.A%i" %(Z,A),
                                      verbose=False,add_info_from_decay_file=True,parallel_run_nb=1,
                                      zip_file_name_decays="../ensdf_data/september_2013/decays.zip",
                                     zip_file_name_levels="../ensdf_data/august_2012/all_levels.zip")
        
    
    print ComputeDecayProductsFromG4Database(Z,A,elevel=0.,database_directory="../g4data/RadioactiveDecay4.0",
                                        photo_evap_database_directory="../g4data/PhotonEvaporation3.0")
    """
"""    
level_energy_vec,level_half_life_vec,level_is_stable_vec,ref_level_vec,level_magnetic_dipole_moment_str,level_Jpi_str=ReadAdoptedLevelInENSDFZipFile(47,129,zip_file_name="../ensdf_data/august_2012/all_levels.zip")
print   level_energy_vec
print level_half_life_vec
print ref_level_vec
"""

    
     