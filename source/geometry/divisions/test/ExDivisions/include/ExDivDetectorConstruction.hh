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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExDivDetectorConstruction_h
#define ExDivDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ExVDivTester.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExDivDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExDivDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    ExDivDetectorConstruction( const G4String& solidTypeStr,
                               const G4String& PVTypeStr,
                               const G4String& PosTypeStr,
                               const std::vector<G4String>& extraPars );
    ~ExDivDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     void ConstructSDandField();
     ExVDivTester* CreateSolidTester( const G4String& stype,
                                      const G4String& thePVTypeStr,
                                      const G4String& thePosTypeStr,
                                      std::vector<G4String>& extraPars );
     PVType    getPVType ( const G4String& pvt );
     PlaceType getPosType( const G4String& pos );
     G4double GetWorldLengthXY()  {return theDivTester->GetWorldLengthXY();}
     G4double GetWorldLengthZ()   {return theDivTester->GetWorldLengthZ();}
     G4double GetWorldGap()       {return theDivTester->GetWorldGap();}

  private:

    G4String theSolidTypeStr;
    G4String thePVTypeStr;
    G4String thePosTypeStr;
    std::vector<G4String> theExtraPars;

    ExVDivTester* theDivTester;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
