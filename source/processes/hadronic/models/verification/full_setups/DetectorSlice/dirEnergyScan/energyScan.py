#!/usr/bin/python

#--------------------------------------------------------------------------
# Last update: 18-Aug-2008
#
# This script requires no input arguments, but it has some parameters that
# need to be specified (see "***LOOKHERE***" below), and should be run as:
#
#         $  python energyScan.py
#
# The script examines all the files of the form:
#
#         output.log-*
#
# in the current directory (so it is a good idea to copy in a single
# directory all the files obtained with the same Geant4 version,
# for all Physics Lists), and produces in output the file:
#
#         energyScan.txt
#
# that contains the following information:
#
#          LHEP    ...  FTFP_BERT
#    1GeV  xxx1_A  ...  xxx1_B
#    1GeV  ...
#    1GeV  yyy1_A  ...  yyy1_B
#    2GeV  xxx2_A  ...  xxx2_B
#    2GeV  ...
#    2GeV  yyy2_A  ...  yyy2_B
#    ...
#   30GeV  xxx30_A ...  xxx30_B
#   30GeV  ...
#   30GeV  yyy30_A ...  yyy30_B
#
# where  xxx1_A  is the mean value of the first observable listObservables[0]
# for the first Physics List listPhysicsLists[0] , as found in the
# corresponding log file:
#   output.log-CALORIMETER-PARTICLE-1GeV-listBeamEnergies[0]-listPhysicsList[0]
# Similarly, xxx1_B  is the mean value of the first observable specified
# in  listObservables[0]  for the last Physics List specified
# in  listPhysicsLists , always for the 1 GeV case.
# yyy1_A  is the mean value of the last observable in listObservable
# for the Physics List  listPhysicsLists[0] , always for the 1 GeV case.
# yyy2_B  is the mean value of the last observables in listObservable
# for the last Physics List  in  listPhysicsLists,
# always for the 1 GeV case.
# Then the same is repeated for all the beam energies, as specified
# in  listBeamEnergies . 
# In other words:
#   - the number of columns of numbers is equal to the number of Physics Lists;
#   - the number of rows of numbers, for a given beam energy, is equal to
#   - the number of observables.
#
# This script uses the Python script  addInfoInTable.py , which in turn
# uses the Python script  printInfoLogfile.py : both of them must be
# present in the same directory as  energyScan.py .
#
# NB) The scripts energyScan.py , addInfoInTable.py , printInfoLogfile.py
#     have been obtained by modifying the Python scripts with the same
#     names in  StatAccepTest/dirEnergyScan/
#
#--------------------------------------------------------------------------

import os
import sys
import string
import math

#***LOOKHERE***

listPhysicsLists = ( \
                     'LHEP',
                     'QGSP',
                     'QGSC',
                     'FTFP',
                     'QGSP_BIC',
                     'QGSP_BERT',
                     'FTF_BIC',
                     'QGS_BIC',
                     'FTFP_BERT',
                   )
                  
#PARTICLE='e-'
PARTICLE='pi-'
#PARTICLE='pi+'
#PARTICLE='p'

#--- 3 possibilities : 1) combinedCMS
#                      2) combinedATLASbarrel 
#                      3) combinedATLASendcap
#                       
CALORIMETER='combinedCMS'
#CALORIMETER='combinedATLASbarrel'
#CALORIMETER='combinedATLASendcap'

listObservables = ( \
###                    'CPU',
###                    'Tracker_E',  
###                    'EM_Evis',
###                    'EM_Etot',
                    'HAD_Evis',
###                    'HAD_Etot',
###                    'Muon_E',
###                    'EM_res',
                    'HAD_res',
###                    'EM_e_Evis',
###                    'EM_mu_Evis',
###                    'EM_pi_Evis',
###                    'EM_k_Evis',
###                    'EM_p_Evis',
###                    'EM_n_Evis',
###                    'HAD_e_Evis',
###                    'HAD_mu_Evis',
###                    'HAD_pi_Evis',
###                    'HAD_k_Evis',
###                    'HAD_p_Evis',
###                    'HAD_n_Evis',
                    )

listBeamEnergies = ( \
                     '1GeV',
                     '2GeV',
                     '3GeV',
                     '4GeV',
                     '5GeV',
                     '6GeV',
                     '7GeV',
                     '8GeV',
                     '9GeV',
                     '10GeV',
                     '11GeV',
                     '12GeV',
                     '13GeV',
                     '14GeV',
                     '15GeV',
                     '16GeV',
                     '17GeV',
                     '18GeV',
                     '19GeV',
                     '20GeV',
                     '21GeV',
                     '22GeV',
                     '23GeV',
                     '24GeV',
                     '25GeV',
                     '26GeV',
                     '27GeV',
                     '28GeV',
                     '29GeV',
                     '30GeV',
                    )

#***endLOOKHERE***


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) != 1 ) :
    print ' Usage:  energyScan.py '
else :
    print ' ------------------------------------------ '
    print ' Parameters chosen: '
    print '    Particle=', PARTICLE
    print '    Calorimeter=', CALORIMETER
    print '    PhysicsLists=', listPhysicsLists
    print '    BeamEnergies=', listBeamEnergies  
    print ' ------------------------------------------ '

    os.system( "rm -f template.txt" )
    theTemplateFile = open( "template.txt", "w" )
    if ( theTemplateFile ) :
        for item in listObservables :
            #print item
            theTemplateFile.write( item + "\n" )
    theTemplateFile.close()

    os.system( "rm -f listLogFiles.txt" )
    os.system( "ls -1 output.log-" + CALORIMETER + "-" + PARTICLE + "* > listLogFiles.txt" )

    stringAllPhysicsLists = ""
    for iPhysicsList in listPhysicsLists :
        stringAllPhysicsLists += iPhysicsList + " "

    os.system( "rm -f energyScan.txt" )
    theResultFile = open( "energyScan.txt", "w" )
    theResultFile.write( "           " + stringAllPhysicsLists + "\n" )

    # The idea is to use the script  addInfoInTable.py  for each
    # beam energy value. In other words, for each beam energy value
    # we find the log-file associated with each Physics List, then
    # we symbolic link this file with the name of the Physics List,
    # and then we call  addInfoInTable.py .
    # Notice that if one of the expected log-files is missing,
    # for any reason, the script  printInfoLogfile.py , which is
    # invoked by  addInfoInTable.py , will fill the missing
    # information with  "0" .

    for iEnergy in listBeamEnergies :
        #print iEnergy
        for iPhysicsList in listPhysicsLists :
            #print iPhysicsList
            os.system( "rm -f " + iPhysicsList )
            theListLogFiles = open( "listLogFiles.txt", "r" )
            if ( theListLogFiles ) :
                for iFile in theListLogFiles :
                    # We have to be careful of not confusing:
                    #   - e.g. FTFP and FTFP_BERT;
                    #          QGSP, QGSP_BERT, and QGSP_BIC;
                    #     when in reality only one should be chosen;
                    #   - e.g. p, pi-, pi+ when in reality only one
                    #     should be chosen
                    #     (remember that we use the following
                    #      inconsistent convention for particles
                    #      names in the output files:
                    #       e.g. output.log-cms-e-1GeV-...   for e- 
                    #            output.log-cms-pi-1GeV-...  for pi-
                    #            output.log-cms-pi+-1GeV-... for pi+
                    #            output.log-cms-p-1GeV-...   for p
                    if ( iFile.rstrip().find( "-" + iEnergy + "-" ) > -1        and 
                         iFile.rstrip().find( "-" + iPhysicsList ) > -1         and 
                         iFile.rstrip().find( "-" + iPhysicsList + "_" ) == -1  and 
                         iFile.rstrip().find( "-" + PARTICLE ) > -1             and 
                         ( iFile.rstrip().find ( "-" + PARTICLE + "-" ) > -1  or
                           PARTICLE.find( "-" ) > -1 )                          and
                         iFile.rstrip().find( "-" + CALORIMETER + "-" ) > -1 ) :
                        #print iFile
                        command = "ln -sfn " + iFile.rstrip() + " " + iPhysicsList
                        #print 'command=', command
                        os.system( command )
            theListLogFiles.close()            

        os.system( "cp template.txt table.txt" )
        command = "python addInfoInTable.py table.txt " + stringAllPhysicsLists
        #print 'command=', command
        os.system( command )

        theTable = open( "table.txt", "r" )
        for line in theTable :
            theResultFile.write( iEnergy + " " + line )
        theTable.close()

        for iPhysicsList in listPhysicsLists :
            os.system( "rm -f " + iPhysicsList )

    theResultFile.close()

    # Clean up temporary files.
    os.system( "rm -f listLogFiles.txt" )
    os.system( "rm -f template.txt" )
    os.system( "rm -f table.txt" )
    os.system( "rm -f old.*" )
    os.system( "rm -f .print*" )

#-------------------------------------------------------------------
