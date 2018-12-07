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
/// \file persistency/P02/include/ExP02DetectorConstruction.hh
/// \brief Definition of the ExP02DetectorConstruction class
//
//
#ifndef ExP02DetectorConstruction_H
#define ExP02DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4Element.hh"

/// Detector construction

class ExP02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    ExP02DetectorConstruction();
    ~ExP02DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

  private:
    
    // Logical volumes
    //
    G4LogicalVolume* fExperimentalHall_log;
    G4LogicalVolume* fTracker_log;
    G4LogicalVolume* fCalorimeterBlock_log;
    G4LogicalVolume* fCalorimeterLayer_log;

    // Physical volumes
    //
    G4VPhysicalVolume* fExperimentalHall_phys;
    G4VPhysicalVolume* fCalorimeterLayer_phys;
    G4VPhysicalVolume* fCalorimeterBlock_phys;
    G4VPhysicalVolume* fTracker_phys;
};

#endif

