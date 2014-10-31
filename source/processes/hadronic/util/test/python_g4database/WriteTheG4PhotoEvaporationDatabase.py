## Python script to write the G4Photoevaporation database 
#author:  L.Desorgher
#
#History: 
#---------
#      26/07/2014     Creation

#import needed modules
import G4Database
import ENSDF
import NuclearWalletCards
import QCalc
import  sys
import Bricc


#Define the location of the different  data needed to build the database
ENSDF.SetZipFilesENSDF("../ensdf_data/september_2013/decays.zip",
                       "../ensdf_data/august_2012/all_levels.zip") 
ENSDF.SetZipFileXUNDL("../ensdf_data/xundl.zip")

NuclearWalletCards.SetNuclearWalletFileName("../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt")
QCalc.SetQCALCDirectory("../qcalc_data/august2012/")

Bricc.SetBriccDirectory("../extra_executable/Bricc_V23")





#Zmin=int(sys.argv[1])
#Zmax=int(sys.argv[2])
#parallel_run_nb=int(sys.argv[3])
Zmin=23
Zmax=23
parallel_run_nb=None
"""
G4Database.WriteG4PhotoEvaporationFile(23,54,"test.txt")
"""

G4Database.WriteAllG4PhotoEvaporationDatabase(new_directory="../final_g4data/test_photo_evap",
                                     Zmin=Zmin,Zmax=Zmax,parallel_run_nb=parallel_run_nb)
