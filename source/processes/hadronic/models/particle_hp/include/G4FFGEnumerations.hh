//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
 * File:   G4FFGEnumerations.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 6, 2011, 9:12 AM
 */

#ifndef G4FFGENUMAERATIONS_HH
#define	G4FFGENUMAERATIONS_HH

#include "G4Types.hh"

/** G4FFGEnumerations is a namespace that contains all the enumerations that
 *  are used in the fission fragment generator code.
 */
namespace G4FFGEnumerations
{
    /** The two types of fission data available. Independent yields are taken
     *  directly from  prompt fission products, while cumulative fission yields are
     *  taken some time after the fission event so that the products have had a
     *  moment to stabilize.
     */
    enum YieldType { INDEPENDENT = 454,
                     CUMULATIVE = 459 };
    /** The first value of YieldType */
    static const G4int YieldTypeFirst = INDEPENDENT;
    /** The last value of YieldType */
    static const G4int YieldTypeLast = CUMULATIVE;

    /** Causes of fission.
     *  Currently only yields for spontaneous and neutron induced fission are
     *  included in the data libraries. Photon and gamma induced fission are
     *  included here only to provide potential forward compatibility.
     *  <b> G4FissionFragmentGenerator::G4SetCause() must be modified if the data
     *  sets are expanded to include photon and gamma induced fission. </b>
     */
    enum FissionCause { SPONTANEOUS,
                        NEUTRON_INDUCED,
                        PROTON_INDUCED,
                        GAMMA_INDUCED };
    /** The first value of FissionCause */
    static const G4int FissionCauseFirst = SPONTANEOUS;
    /** The last value of FissionCause.
      * This is set to NEUTRON_INDUCED becuase neither PROTON_INDUCED 
      * nor GAMMA_INDUCED are currently supporded.
      */
    static const G4int FissionCauseLast = NEUTRON_INDUCED;

    /** The possible fission sampling methods */
    enum FissionSamplingScheme{ NORMAL,
                                LIGHT_FRAGMENT };
    /** The first value of FissionSamplingScheme */
    static const G4int FissionSamplingSchemeFirst = NORMAL;
    /** The last value of FissionSamplingScheme */
    static const G4int FissionSamplingSchemeLast = LIGHT_FRAGMENT;

    /** Truncate the Gaussian distribution at 0 (\p POSITIVE) or sample all values
     *  (\p ALL)
     */
    enum GaussianRange { POSITIVE,
                         ALL };

    /** Sample a discretized Gaussian distribution (\p INT) or continuous (\p DOUBLE) */
    enum GaussianReturnType { INT,
                              DOUBLE };

    /** ENDF format provides for 3 isomers - 1 ground state and 2 meta states */
    enum MetaState { GROUND_STATE,
                     META_1,
                     META_2 };
    /** The first value of MetaState */
    static const G4int MetaStateFirst = GROUND_STATE;
    /** The last value of MetaState */
    static const G4int MetaStateLast = META_2;

    /** These are the source shapes available */
    enum SourceType { RECTANGLE,
                      CYLINDER,
                      SPHERE };
    /** The first value of SourceType */
    static const G4int SourceTypeFirst = RECTANGLE;
    /** The last value of SourceType */
    static const G4int SourceTypeLast = SPHERE;

    /** These are the verbosity levels */
    enum Verbosity { SILENT = 0x0,
                     UPDATES = 0x1,
                     DAUGHTER_INFO = 0x2,
                     NEUTRON_INFO = 0x4,
                     GAMMA_INFO = 0x8,
                     ALPHA_INFO = 0x10,
                     MOMENTUM_INFO = 0x20,
                     ENERGY_INFO = 0x40,
                     EXTRAPOLATION_INTERPOLATION_INFO = 0x80,
                     WARNING = 0x100,
                     DEBUG = 0x200,
                     PRINT_ALL = 0x3FF,
                     REPRESS_FUNCTION_ENTER_LEAVE_MESSAGES = 0x1000,
                     REPRESS_RANDOM_SAMPLING_MESSAGES = 0x2000,
                     REPRESS_RECURSIVE_DEBUG_MESSAGES = 0x4000,
                     REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES = 0x8000 };
    /** This is the maximum value of Verbosity */
    static const G4int VerbosityLast = (REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES << 1) - 1;
}

#endif	/* G4FFGENUMAERATIONS_HH */

