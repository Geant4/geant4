#!/usr/bin/python

#--------------------------------------------------------------------------
# Last update: 01-Aug-2008
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

#--- Only 4 possibilities : 1) cms :   Cu-Sci  10-lambda simplified CMS HCAL calorimeter
#                           2) atlas : Cu-LAr  10-lambda simplified ATLAS HEC calorimeter
#                           3) pbwo4 : PbWO4   10-lambda simplified calorimeter
#                           4) tile :  Fe-Sci  10-lambda simplified ATLAS TILECAL calo.
CALORIMETER='cms'
###CALORIMETER='atlas'
###CALORIMETER='pbwo4'
###CALORIMETER='tile'

listObservables = ( \
###                    'CPU',
                    'Evis',
###                    'Etot',
###                    'fL1',
###                    'fL2',
###                    'fL3',
###                    'fL4',
###                    'fR1',
###                    'fR2',
###                    'fR3',
                    'res',
###                    'e_Evis',
###                    'e_Tot',
###                    'e_fL1',
###                    'e_fL2',
###                    'e_fL3',
###                    'e_fL4',
###                    'e_fR1',
###                    'e_fR2',
###                    'e_fR3',
###                    'mu_Evis',
###                    'mu_Etot',
###                    'pi_Evis',
###                    'pi_Etot',
###                    'pi_fL1',
###                    'pi_fL2',
###                    'pi_fL3',
###                    'pi_fL4',
###                    'pi_fR1',
###                    'pi_fR2',
###                    'pi_fR3',
###                    'k_Evis',
###                    'k_Etot',
###                    'p_Evis',
###                    'p_Etot',
###                    'p_fL1',
###                    'p_fL2',
###                    'p_fL3',
###                    'p_fL4',
###                    'p_fR1',
###                    'p_fR2',
###                    'p_fR3',
###                    'n_Evis',
###                    'n_Etot',
###                    'n_fL1',
###                    'n_fL2',
###                    'n_fL3',
###                    'n_fL4',
###                    'n_fR1',
###                    'n_fR2',
###                    'n_fR3',
###                    '#steps_tot',
###                    '#steps_+',
###                    '#steps_0',
###                    '#steps_-',
###                    '#steps_nuclei',
###                    '#steps_unknown',
###                    '#steps_EM',
###                    '#steps_EWK',
###                    '#steps_had',
###                    '#steps_mes',
###                    '#steps_bar',
###                    '#steps_lmes',
###                    '#steps_lbar',
###                    '#steps_smes',
###                    '#steps_sbar',
###                    '#steps_hmes',
###                    '#steps_hbar',
###                    '#steps_e-',
###                    '#steps_gam',
###                    '#steps_e+',
###                    '#steps_mu-',
###                    '#steps_mu+',
###                    '#steps_tau-',
###                    '#steps_tau+',
###                    '#steps_nu',
###                    '#steps_pi+',
###                    '#steps_pi0',
###                    '#steps_pi-',
###                    '#steps_k+',
###                    '#steps_k0',
###                    '#steps_k-',
###                    '#steps_p',
###                    '#steps_pbar',
###                    '#steps_n',
###                    '#steps_nbar',
###                    '#tot',
###                    '#+',
###                    '#0',
###                    '#-',
###                    '#nuclei',
###                    '#unknown',
###                    '#EM',
###                    '#EWK',
###                    '#had',
###                    '#mes',
###                    '#bar',
###                    '#lmes',
###                    '#lbar',
###                    '#smes',
###                    '#sbar',
###                    '#hmes',
###                    '#hbar',
###                    '#e-',
###                    '#gam',
###                    '#e+',
###                    '#mu-',
###                    '#mu+',
###                    '#tau-',
###                    '#tau+',
###                    '#nu',
###                    '#pi+',
###                    '#pi0',
###                    '#pi-',
###                    '#k+',
###                    '#k0',
###                    '#k-',
###                    '#p',
###                    '#pbar',
###                    '#n',
###                    '#nbar',
###                    'L_e',
###                    'L_pi',
###                    'L_p',
###                    'L_gamma',
###                    'L_n',
###                    'exit_kin',
###                    'exit_fgam',
###                    'exit_fn',
###                    'exit_fnu',
###                    'exit_fmu',
###                    'exit_fe',
###                    'exit_fother',
###                    'exit_#',
###                    'exit_#gam',
###                    'exit_#n',
###                    'exit_#nu',
###                    'exit_#mu',
###                    'exit_#e',
###                    'exit_#other"',
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
            os.system( "rm -f iPhysicsList" )
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
