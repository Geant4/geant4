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
// $Id: HadrontherapyDetectorConstruction.hh,v 3.0, September 2004
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G. Candiano, G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapyModulator_H
#define HadrontherapyModulator_H 1

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

 private:
  const G4int numberOfModules;
  G4VPhysicalVolume* physiMotherMod;
  G4Tubs** solidMod;
  G4LogicalVolume** logicMod;
  G4VPhysicalVolume** physiMod; 
  G4String file;
  G4DataVector* startAngle;
  G4DataVector* spanningAngle;
  G4DataVector* Ztranslation;
  G4RotationMatrix* rotationMatrix;
   
  void ReadFile(G4String);
  void ReadData(G4String);
};
#endif
