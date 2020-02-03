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
/// \file Par02DetectorConstruction.hh
/// \brief Definition of the Par02DetectorConstruction class

#ifndef PAR02_DETECTOR_CONSTRUCTION_H
#define PAR02_DETECTOR_CONSTRUCTION_H

#include "Par02FastSimModelTracker.hh"
#include "Par02FastSimModelEMCal.hh"
#include "Par02FastSimModelHCal.hh"

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"
#include "G4GlobalMagFieldMessenger.hh"

/// Construction of detector geometry.
///
/// A mandatory initialization class of the detector setup. 
/// Detector construction allows to use the geometry read from a GDML file. 
/// Based on G4 examples/persistency/gdml/G01/include/G01DetectorConstruction.hh .
/// @author Anna Zaborowska

class Par02DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    
    /// A default constructor.
    Par02DetectorConstruction();

    virtual ~Par02DetectorConstruction();
    
    /// A method invoked by G4RunManager::Initialize()
    /// @return A pointer to the world volume.
    virtual G4VPhysicalVolume* Construct();

    /// A method invoked by G4RunManager::Initialize() to construct thread local objects
    virtual void ConstructSDandField();

    /// A vector of the tracking detector regions
    std::vector< G4Region* > fTrackerList;

    /// A vector of the the electromagnetic calorimeter regions
    std::vector< G4Region* > fECalList;

    /// A vector of the the hadronic calorimeter regions
    std::vector< G4Region* > fHCalList;

    /// A vector of the muon detector regions
    std::vector< G4Region* > fMuonList;

    /// Messenger of the magnetic field
    G4GlobalMagFieldMessenger* fMagFieldMessenger;
};

#endif

