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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyModulator_H
#define HadrontherapyModulator_H 1

#include "globals.hh"
#include <fstream>
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "HadrontherapyModulatorMessenger.hh"

// using namespace std;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;


class HadrontherapyModulator
{
public:

  HadrontherapyModulator();
  ~HadrontherapyModulator();

  void BuildModulator(G4VPhysicalVolume*);  
  void SetModulatorAngle(G4double);
  void SetModulatorMaterial(G4String);
  void SetModulatorPosition(G4ThreeVector);
  void SetModulatorInnerRadius(G4double);
  void SetModulatorOuterRadius(G4double);
  void  ModulatorDefaultProperties();
  void ModulatorPropertiesFromFile(G4String);
  void GetDataFromFile(G4String value);
  void  GetStepInformation();
  void BuildSteps();

private:
  std::ifstream File;
   
   G4LogicalVolume * logicMotherMod ;
   G4VPhysicalVolume* physiMotherMod;
   
   G4Material* Mod0Mater;
  G4Material* ModMater;

  G4Tubs*           solidMod1;    
  G4LogicalVolume*  logicMod1;   
  G4VPhysicalVolume* physiMod1;
  
  G4Tubs*            solidMod2;   
  G4LogicalVolume*   logicMod2;   
  G4VPhysicalVolume* physiMod2;
  
  G4Tubs*            solidMod3;   
  G4LogicalVolume*   logicMod3;   
  G4VPhysicalVolume* physiMod3;
  
  G4Tubs*            solidMod4;  
  G4LogicalVolume*   logicMod4;  
  G4VPhysicalVolume* physiMod4;
  
  G4double pi;
  G4int StepNumbers;
  G4double* Weight;
  G4double* StepThickness;
  G4double* StartingAngle;
  G4double* SpanningAngle;
  G4ThreeVector* PositionMod;
  G4Tubs** solidMod;
  G4LogicalVolume**	logicMod;
  G4VPhysicalVolume**	physiMod;

  G4RotationMatrix* rm;
  
  G4String FileName;
  HadrontherapyModulatorMessenger* ModulatorMessenger;
  G4double innerRadiusOfTheTube ;
  G4double outerRadiusOfTheTube ;
 

   
};
#endif
