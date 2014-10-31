## Python script to write the G4Radecay database 
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

#Zmin=int(sys.argv[1])
#Zmax=int(sys.argv[2])


#Define the location of the different  data needed to build the database
ENSDF.SetZipFilesENSDF("../ensdf_data/september_2013/decays.zip","../ensdf_data/august_2012/all_levels.zip") 
ENSDF.SetZipFileXUNDL("../ensdf_data/xundl.zip")

NuclearWalletCards.SetNuclearWalletFileName("../nudat2_data/nuclear_wallet/august2012/all_nuclei.txt")
QCalc.SetQCALCDirectory("../qcalc_data/august2012/")




Zmin=12
Zmax=14
G4Database.WriteAllG4RadDecayDatabase(new_directory="../final_g4data/test",Zmin=Zmin,Zmax=Zmax)