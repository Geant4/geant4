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
//
/// \file Par02DetectorParametrisation.hh
/// \brief Definition of the Par02DetectorParameterisation class

#ifndef PAR02_DETECTOR_PARAMETRSIATION_H
#define PAR02_DETECTOR_PARAMETRSIATION_H

#include "globals.hh"

/// Definition of detector resolution and efficiency.
///
/// A simple class used to provide the detector resolution and efficiency
/// (dependent on the detector, parametrisation type and particle momentum).
/// @author Anna Zaborowska

class Par02DetectorParametrisation {
  public:
    
    /// A default constructor.
    Par02DetectorParametrisation();

    ~Par02DetectorParametrisation();
    
    /// A parametrisation type (CMS, ATLAS, ALEPH).
    enum Parametrisation { eCMS, eATLAS, eALEPH };
    
    /// A detector type (tracking detector, electromagnetic calorimeter,
    ///                  hadronic calorimeter).
    enum Detector { eTRACKER, eEMCAL, eHCAL };
    
    /// Gets the resolution of a detector for a given particle.
    /// @param aDetector A detector type.
    /// @param aParametrisation A parametrisation type.
    /// @param aMomentum A particle momentum.
    G4double GetResolution( Detector aDetector, Parametrisation aParametrisation,
                            G4double aMomentum );
    
    /// Gets the efficiency of a detector for a given particle.
    /// @param aDetector A detector type.
    /// @param aParametrisation A parametrisation type.
    /// @param aMomentum A particle momentum.
    G4double GetEfficiency( Detector aDetector, Parametrisation aParametrisation,
                            G4double aMomentum );
};

#endif

