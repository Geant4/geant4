#!/usr/bin/python

#--------------------------------------------------------------------------
# Last update: 05-May-2009
#
# This script requires at least 2 input arguments, and should be
# run as:
#
#         $  python addInfoInTable.py table file1 [file2] ... [fileN]
#
# where the  "table"  is a text (ascii) file with a table, where
# at least the first column, with the names of the observables,
# should be present, and the second argument, "file1" should be
# a log-files produced by running one of the StatAccepTest
# simulations, of Fluka simulations.
# The remain arguments,  file2 ... fileN  are optionals,
# and they are also meant to be log-files of StatAccepTest simulations,
# or Fluka simulations.
# The table and files must be in the same directory as the script
# (use symbolic link if they are located in other directories).
#
# This script produces in output a new table, with the same name
# as the one given in input (whereas the original one is backup
# with name "old.0." in front of the original name; all other
# eventual intermediate tables are saved with names "old.1.",
# "old.2.", etc., in front of the original name: this is useful
# for debugging), in which new columns have been added with the
# corresponding information extract from the input files:
#    file1  first,  file2  second,  etc...
# using the script  printInfoLogFile.py  (which should be
# located in the same directory as addInfoInTable.py ).
# (For debugging, the output of  printInfoLogFile.py , when run
#  with each of the input files, can be found in:
#  .printInfoLogfile.out-file1 , .printInfoLogfile.out-file2, etc. )
#
# Notice that the rows corresponding to observables which are
# not recognized are reproduced in the new table as they are.
# If one of the files specified in the argument is not a log-file
# obtained by running one of the StatAccepTest simulations, or
# Fluka simulatiosn, nothing will happen, i.e. the table will not
# be updated.
#
# Look ***LOOKHERE*** to map the short names of the observables,
# used in the first column of the table, with the actual names
# of the observables as they appear in "file".
# However, if you have to build a new table, these are the names
# of the observables that can be used directly as the first column
# of the table, without requiring any change in the script below
# (of course, you can use a subset of them):
#
#        CPU
#        Evis
#        Etot
#        fL1
#        fL2
#        fL3
#        fL4
#        fR1
#        fR2
#        fR3
#        res
#        e_Evis
#        e_Tot
#        e_fL1
#        e_fL2
#        e_fL3
#        e_fL4
#        e_fR1
#        e_fR2
#        e_fR3
#        mu_Evis
#        mu_Etot
#        pi_Evis
#        pi_Etot
#        pi_fL1
#        pi_fL2
#        pi_fL3
#        pi_fL4
#        pi_fR1
#        pi_fR2
#        pi_fR3
#        k_Evis
#        k_Etot
#        p_Evis
#        p_Etot
#        p_fL1
#        p_fL2
#        p_fL3
#        p_fL4
#        p_fR1
#        p_fR2
#        p_fR3
#        n_Evis
#        n_Etot
#        n_fL1
#        n_fL2
#        n_fL3
#        n_fL4
#        n_fR1
#        n_fR2
#        n_fR3
#        #steps_tot
#        #steps_+
#        #steps_0
#        #steps_-
#        #steps_nuclei
#        #steps_unknown
#        #steps_EM
#        #steps_EWK
#        #steps_had
#        #steps_mes
#        #steps_bar
#        #steps_lmes
#        #steps_lbar
#        #steps_smes
#        #steps_sbar
#        #steps_hmes
#        #steps_hbar
#        #steps_e-
#        #steps_gam
#        #steps_e+
#        #steps_mu-
#        #steps_mu+
#        #steps_tau-
#        #steps_tau+
#        #steps_nu
#        #steps_pi+
#        #steps_pi0
#        #steps_pi-
#        #steps_k+
#        #steps_k0
#        #steps_k-
#        #steps_p
#        #steps_pbar
#        #steps_n
#        #steps_nbar
#        #tot
#        #+
#        #0
#        #-
#        #nuclei
#        #unknown
#        #EM
#        #EWK
#        #had
#        #mes
#        #bar
#        #lmes
#        #lbar
#        #smes
#        #sbar
#        #hmes
#        #hbar
#        #e-
#        #gam
#        #e+
#        #mu-
#        #mu+
#        #tau-
#        #tau+
#        #nu
#        #pi+
#        #pi0
#        #pi-
#        #k+
#        #k0
#        #k-
#        #p
#        #pbar
#        #n
#        #nbar
#        L_e
#        L_pi
#        L_p
#        L_gamma
#        L_n
#        exit_kin
#        exit_fgam
#        exit_fn
#        exit_fnu
#        exit_fmu
#        exit_fe
#        exit_fother
#        exit_#
#        exit_#gam
#        exit_#n
#        exit_#nu
#        exit_#mu
#        exit_#e
#        exit_#other
#
# NB) In the case of Fluka simulations some of the above observables
#     (CPU, #steps, most of the # of particles, L_*, and some of
#      exit informations) are not defined. However, you can run this
#     script normally because it will be printed "0" for those
#     observables that are not defined in Fluka.
#
#--------------------------------------------------------------------------

import os
import sys
import string
import math


#===============================================
#================= FUNCTIONS =================== 
#===============================================


def updateLine( infoList, tableLine ) :
    # This function received in input two arguments:
    # the first is a list whose elements are the lines of the
    # file with the information we want to add;
    # the second one is a line of the table which we want to
    # update with the information provided in the file.
    # This function returns the new, updated line.

    #***LOOKHERE*** : here is the map between the short names of the
    #                 observables (key of the dictionary), and the
    #                 complete name as it appears in "infoFile"
    #                 (value of the dictionary).
    mapNameObservables = {
        #
        'CPU':'cpuTime',
        'Evis':'visEnergy',
        'Etot':'totEnergy',
        'fL1':'fL1',
        'fL2':'fL2',
        'fL3':'fL3',
        'fL4':'fL4',
        'fR1':'fR1',
        'fR2':'fR2',
        'fR3':'fR3',
        'res':'resolution',
        'e_Evis':'electronEvis',
        'e_Tot':'electronEtot',
        'e_fL1':'electronfL1',
        'e_fL2':'electronfL2',
        'e_fL3':'electronfL3',
        'e_fL4':'electronfL4',
        'e_fR1':'electronfR1',
        'e_fR2':'electronfR2',
        'e_fR3':'electronfR3',
        'mu_Evis':'muonEvis',
        'mu_Etot':'muonEtot',
        'pi_Evis':'pionEvis',
        'pi_Etot':'pionEtot',
        'pi_fL1':'pionfL1',
        'pi_fL2':'pionfL2',
        'pi_fL3':'pionfL3',
        'pi_fL4':'pionfL4',
        'pi_fR1':'pionfR1',
        'pi_fR2':'pionfR2',
        'pi_fR3':'pionfR3',
        'k_Evis':'kaonEvis',
        'k_Etot':'kaonEtot',
        'p_Evis':'protonEvis',
        'p_Etot':'protonEtot',
        'p_fL1':'protonfL1',
        'p_fL2':'protonfL2',
        'p_fL3':'protonfL3',
        'p_fL4':'protonfL4',
        'p_fR1':'protonfR1',
        'p_fR2':'protonfR2',
        'p_fR3':'protonfR3',
        'n_Evis':'nucleusEvis',
        'n_Etot':'nucleusEtot',
        'n_fL1':'nucleusfL1',
        'n_fL2':'nucleusfL2',
        'n_fL3':'nucleusfL3',
        'n_fL4':'nucleusfL4',
        'n_fR1':'nucleusfR1',
        'n_fR2':'nucleusfR2',
        'n_fR3':'nucleusfR3',
        '#steps_tot':'# total  steps',
        '#steps_+':'# positives  steps',
        '#steps_0':'# neutrals  steps',
        '#steps_-':'# negatives  steps',
        '#steps_nuclei':'# nuclei  steps',
        '#steps_unknown':'# particles with Unrecognized PDG code  steps',
        '#steps_EM':'# electromagnetic (e+ , e- , gammas)  steps',
        '#steps_EWK':'# electroweak (mu+, mu-, tau+, tau-, neutrinos)  steps',
        '#steps_had':'# hadrons  steps',
        '#steps_mes':'# mesons  steps',
        '#steps_bar':'# baryons  steps',
        '#steps_lmes':'# light mesons (u/ubar/d/dbar)  steps',
        '#steps_lbar':'# light baryons (u/ubar/d/dbar)  steps',
        '#steps_smes':'# strange (s/sbar) mesons  steps',
        '#steps_sbar':'# strange (s/sbar) baryons  steps',
        '#steps_hmes':'# heavy (c/cbar or b/bbar) mesons  steps',
        '#steps_hbar':'# heavy (c/cbar or b/bbar) baryons  steps',
        '#steps_e-':'# electrons  steps',
        '#steps_gam':'# gammas  steps',
        '#steps_e+':'# positrons  steps',
        '#steps_mu-':'# mu-  steps',
        '#steps_mu+':'# mu+  steps',
        '#steps_tau-':'# tau-  steps',
        '#steps_tau+':'# tau+  steps',
        '#steps_nu':'# neutrinos  steps',
        '#steps_pi+':'# pi+  steps',
        '#steps_pi0':'# pi0  steps',
        '#steps_pi-':'# pi-  steps',
        '#steps_k+':'# K+  steps',
        '#steps_k0':'# K-neutral (K0/K0bar or K0_S/K0_L)  steps',
        '#steps_k-':'# K-  steps',
        '#steps_p':'# protons  steps',
        '#steps_pbar':'# anti-protons  steps',
        '#steps_n':'# neutrons  steps',
        '#steps_nbar':'# anti-neutrons  steps',
        '#tot':'# total  tracks',
        '#+':'# positives  tracks',
        '#0':'# neutrals  tracks',
        '#-':'# negatives  tracks',
        '#nuclei':'# nuclei  tracks',
        '#unknown':'# particles with Unrecognized PDG code  tracks',
        '#EM':'# electromagnetic (e+ , e- , gammas)  tracks',
        '#EWK':'# electroweak (mu+, mu-, tau+, tau-, neutrinos)  tracks',
        '#had':'# hadrons  tracks',
        '#mes':'# mesons  tracks',
        '#bar':'# baryons  tracks',
        '#lmes':'# light mesons (u/ubar/d/dbar)  tracks',
        '#lbar':'# light baryons (u/ubar/d/dbar)  tracks',
        '#smes':'# strange (s/sbar) mesons  tracks',
        '#sbar':'# strange (s/sbar) baryons  tracks',
        '#hmes':'# heavy (c/cbar or b/bbar) mesons  tracks',
        '#hbar':'# heavy (c/cbar or b/bbar) baryons  tracks',
        '#e-':'# electrons  tracks',
        '#gam':'# gammas  tracks',
        '#e+':'# positrons  tracks',
        '#mu-':'# mu-  tracks',
        '#mu+':'# mu+  tracks',
        '#tau-':'# tau-  tracks',
        '#tau+':'# tau+  tracks',
        '#nu':'# neutrinos  tracks',
        '#pi+':'# pi+  tracks',
        '#pi0':'# pi0  tracks',
        '#pi-':'# pi-  tracks',
        '#pi':' # pi  tracks',
        'pi0/pi':'pi0/pi',
        '#k+':'# K+  tracks',
        '#k0':'# K-neutral (K0/K0bar or K0_S/K0_L)  tracks',
        '#k-':'# K-  tracks',
        '#p':'# protons  tracks',
        '#pbar':'# anti-protons  tracks',
        '#n':'# neutrons  tracks',
        '#nbar':'# anti-neutrons  tracks',
        'L_e':'electronLength',
        'L_pi':'pionLength',
        'L_p':'protonLength',
        'L_gamma':'gammaLength',
        'L_n':'neutronLength',
        'exit_kin':'exitKin',
        'exit_fgam':'exitFracGammas',
        'exit_fn':'exitFracNeutrons',
        'exit_fnu':'exitFracNeutrinos',
        'exit_fmu':'exitFracMuons',
        'exit_fe':'exitFracElectrons',
        'exit_fother':'exitFracOthers',
        'exit_#':'exitNum',
        'exit_#gam':'exitNumGammas',
        'exit_#n':'exitNumNeutrons',
        'exit_#nu':'exitNumNeutrinos',
        'exit_#mu':'exitNumMuons',
        'exit_#e':'exitNumElectrons',
        'exit_#other':'exitNumOthers',
        #
        }

    #print ' tableLine=', tableLine

    newTableLine = tableLine

    # The case of "pi0/pi" or "#pi" needs a special treatment,
    # because the number of pi+, pi-, and pi0 are all needed.
    if ( tableLine.find( "pi0/pi" ) > -1  or
         ( tableLine.find( "#pi" ) > -1    and
           tableLine.find( "#pi+" ) == -1  and
           tableLine.find( "#pi0" ) == -1  and
           tableLine.find( "#pi-" ) == -1 ) ) :
        #print 'tableLine=', tableLine
        numPionPlus = 0.0
        numPionMinus = 0.0
        numPionZero = 0.0
        for infoLine in infoList :
            #print 'infoLine=', infoLine
            if ( infoLine.find( "# pi+  tracks" ) > -1 ) :
                numPionPlus = float( infoLine.split( "=" )[1] )
            elif ( infoLine.find( "# pi0  tracks" ) > -1 ) :
                numPionZero = float( infoLine.split( "=" )[1] )
            elif ( infoLine.find( "# pi-  tracks" ) > -1 ) :
                numPionMinus = float( infoLine.split( "=" )[1] )
        #print ' numPionPlus = ', numPionPlus
        #print ' numPionZero = ', numPionZero
        #print ' numPionMinus = ', numPionMinus
        sum = numPionPlus + numPionZero + numPionMinus
        #print ' sum=', sum
        if ( tableLine.find( "pi0/pi" ) > -1  and  sum > 0.0 ) :
            ratio = 100.0 * numPionZero / sum 
            ratioStringFormatted = "%.0f" % ratio
            newTableLine += "   " + ratioStringFormatted + "%"
        elif ( sum > 0.0 ) :
            sumStringFormatted = "%.0f" % sum
            newTableLine += "   " + sumStringFormatted

    else :
        #print ' tableLine.split()=', tableLine.split()
        if ( len( tableLine.split() ) > 0 ) :
            shortName = tableLine.split()[0]
            #print ' shortName=', shortName
            if ( mapNameObservables.has_key( shortName ) ) :
                longName = mapNameObservables[ shortName ]
                #print ' longName=', longName
                for infoLine in infoList :
                    #print 'infoLine=', infoLine
                    if ( infoLine.find( longName ) > -1 ) :
                        value = float( infoLine.split( "=" )[1] )
                        #print ' value=', value

                        if ( longName.find( "cpuTime" ) > -1 or
                             longName.find( "tracks" ) > -1  or
                             longName.find( "steps" ) > -1   or
                             longName.find( "exit" ) > -1    or
                             longName.find( "neutronLength" ) > -1 ) :
                            valueStringFormatted = "%.0f" % value
                            newTableLine += "   " + valueStringFormatted
                        else :
                            newTableLine += "   " + str( value )
                            
                        if ( longName == "resolution"  or
                             longName == "exitFracNeutrons"  or
                             longName.find( "fL" ) > -1  or
                             longName.find( "fR" ) > -1  or
                             longName.find( "nEvis" ) > -1  or
                             longName.find( "nEtot" ) > -1  or
                             longName.find( "sEvis" ) > -1  or
                             longName.find( "sEtot" ) > -1  ) :
                            newTableLine += "%"

                        break

    #print ' newTableLine=', newTableLine
    
    return newTableLine


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) < 3 ) :
    print " Usage:  addInfoInTable.py table file1 [file2] ... [fileN]"
else :

    # Loop over the input files: file1 , file2 ... fileN
    for iFile in xrange( len( sys.argv ) - 2 ) :
        print ' Considering file : ', sys.argv[ iFile+2 ]

        # Use the script  printInfoLogFile.py .
        os.system( "python printInfoLogfile.py " + sys.argv[ iFile+2 ] + \
                   " > .printInfoLogfile.out-" + sys.argv[ iFile+2 ] + " 2>&1 " )
    
        # The file  .printInfoLogfile.out-...  contains the
        # information we need in order to update the table.
        fileInformation = open( ".printInfoLogfile.out-" + sys.argv[ iFile+2 ], "r" )

        # Put the content of the file inside a list, in order to
        # be able to scan it more times.
        listInformation = []
        for infoLine in fileInformation :
            listInformation.append( infoLine.rstrip() )
            #print 'infoLine=', infoLine,

        print '  -> number of useful information = ', len( listInformation )
        #for lineInformation in listInformation :
        #    print 'lineInformation=', lineInformation
        
        fileInformation.close()

        if ( len( listInformation ) > 0 ) :
                
            table = open( sys.argv[1], "r" )
            if ( table ) :
                    
                # Create a new table.
                newtable = open( "new." + sys.argv[1], "w" )
                if ( newtable ) :
                    for line in table :
                        #print line,
                        newline = updateLine( listInformation, line.rstrip() )
                        for item in newline.split() :
                            mystring = item.ljust( 12 )[:10]
                            newtable.write( mystring )
                        newtable.write( "\n" )

                newtable.close()
            table.close()

            # Rename the new table as the original one, but store 
            # also the previous table (for debugging).
            os.system( "mv " + sys.argv[1] + " old." + str( iFile ) + "." +
                       sys.argv[1] + 
                       " ; mv new." + sys.argv[1] + " " + sys.argv[1] )

#-------------------------------------------------------------------
