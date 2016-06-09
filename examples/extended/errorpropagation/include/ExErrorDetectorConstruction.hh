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
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
//
// Class Description:
//
//  Creates a simplified typical HEP detector: 
//    An air beamline ( BEAM )
//    An air central detector ( CDET )
//    A copper calorimeter, divided in four ( ECAL )
//    An aluminium calorimeter, divided in ten ( HCAL )
//    An air muon detector ( MUON )
//
// History:
// - Created:   P. Arce   May 2007
//

#ifndef ExErrorDetectorConstruction_hh
#define ExErrorDetectorConstruction_hh 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ExErrorMagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExErrorDetectorMessenger;
class G4UserLimits;

//------------------------------------------------------------------------
class ExErrorDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  ExErrorDetectorConstruction();
  ~ExErrorDetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  
  void SetMagField(G4double);

private:
  G4double xBEAM;
  G4double xCDET;
  G4double xECAL;
  G4double xSOLN;
  G4double xHCAL;
  G4double xMUON;
  G4double ndivECAL;
  G4double ndivHCAL;
  G4double yzLength;
  G4double xHalfWorldLength;


  G4UserLimits*         userLimits;
  
  ExErrorMagneticField* fpMagField;   // pointer to the magnetic field 
  
  ExErrorDetectorMessenger* detectorMessenger;  // pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
