//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExDivDetectorConstruction.hh,v 1.2 2004-05-13 14:57:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
