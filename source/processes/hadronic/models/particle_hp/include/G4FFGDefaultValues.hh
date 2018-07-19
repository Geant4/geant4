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
 * File:   G4FFGDefaultValues.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on August 10, 2012, 17:03
 */

#ifndef G4FFGDEFAULTVALUES_HH
#define G4FFGDEFAULTVALUES_HH

//#include "G4SystemOfUnits.hh" bad, use instead:
#include <CLHEP/Units/SystemOfUnits.h>
//  and use constants with namespace, see below e.g. eV becomes CLHEP::eV
#include "G4ThreeVector.hh"

#include "G4FFGEnumerations.hh"

/** G4FFGDefaultValues is a one-stop shop for storing the default values to
 *  variables that configure how the fission fragment generator code is
 *  initialized.
 */
namespace G4FFGDefaultValues
{
// Global
    /** The energy of thermal neutrons */
    static const G4double ThermalNeutronEnergy =                            0.0253 * CLHEP::eV;

// Verbosity
#ifdef G4DEBUG_VERBOSE
    /** Verbosity for the entire package */
    static const G4int Verbosity =                                          G4FFGEnumerations::PRINT_ALL;// | G4FFGEnumerations::REPRESS_FUNCTION_ENTER_LEAVE_MESSAGES;
#else /* G4DEBUG_VERBOSE */
    /** Verbosity for the entire package */
    static const G4FFGEnumerations::Verbosity Verbosity =                   G4FFGEnumerations::SILENT;
#endif /* G4DEBUG_VERBOSE */

// Fission Parameters
    /** Default Isotope */
    static const G4int Isotope =                                            92238;
    /** Default meta state */
    static const G4FFGEnumerations::MetaState MetaState =                   G4FFGEnumerations::GROUND_STATE;
    /** Default fission cause */
    static const G4FFGEnumerations::FissionCause FissionCause =             G4FFGEnumerations::SPONTANEOUS;
    /** Default incident energy */
    static const G4double IncidentEnergy =                                  ThermalNeutronEnergy / CLHEP::eV;
    /** Default incident energy unit */
    static const char IncidentEnergyUnit[] =                                "eV";
    /** Default yield type */
    static const G4FFGEnumerations::YieldType YieldType =                   G4FFGEnumerations::INDEPENDENT;
    /** Default sampling scheme */
    static const G4FFGEnumerations::FissionSamplingScheme SamplingScheme =  G4FFGEnumerations::NORMAL;
    /** Default probabilility of a ternary fission */
    static const G4double TernaryProbability =                              0;
    /** Default alpha production in a ternary fission */
    static const G4double AlphaProduction =                                 0;
    /** Default mean gamma energy for gamma sampling */
    static const G4double MeanGammaEnergy =                                 800 * CLHEP::keV;

// Event Parameters
    /** Default event time */
    static const G4double EventTime =                                       0;
    /** Default event time units */
    static const char EventTimeUnit[] =                                     "ns";

// Source Description
    /** Default source center */
    static const G4ThreeVector SourceCenter(0, 0, 0);
    /** Default source depth */
    static const G4double SourceDepth =                                     1;
    /** Default source rectangle Height */
    static const G4double SourceHeight =                                    1;
    /** Default source radius */
    static const G4double SourceRadius =                                    1;
    /** Default source rectangle Width */
    static const G4double SourceWidth =                                     1;
    /** Default event time units */
    static const char SourceDimensionUnit[] =                               "cm";
    /** Default source type */
    static const G4FFGEnumerations::SourceType SourceType =                 G4FFGEnumerations::SPHERE;

// Messenger
    /** Default command directory */
    static const char UICommandDirectory[] =                                "/process/hadronic/ffgupga";

// Data
    /** ENDF data tape location, reference against \p G4HPNEUTRONDATA */
    static const char ENDFFissionDataLocation[] =                           "/Fission/FF/";
}

#endif /** G4FFGDEFAULTVALUES_HH */


