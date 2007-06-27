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
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "globals.hh"

class DetectorMessenger;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction(RunAction*, PrimaryGeneratorAction*);
  virtual ~DetectorConstruction();

public:
  
  void SetWorldMaterial(const G4String&);
  void SetAbsorberMaterial (const G4String&);
  void SetAbsorberWidth (G4double);

  void SetBeamMomentum(G4double val)       {gen->SetMomentum(val);
                                            run->SetProjectileMomentum(val);}
  void SetBeamMomentumSpread(G4double val) {gen->SetMomentumSpread(val);
                                            run->SetMomentumSpread(val);}
  void SetFileName(const G4String& nam)    {run->SetFileName(nam);}

  G4VPhysicalVolume* GetPhysWorld()        {return physWorld;}
  G4VPhysicalVolume* GetPhysAbsorber()     {return physAbs;}

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
     
private:

  G4Material*        matWorld;
  G4Material*        matAbsorber;

  G4double           width;

  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physAbs;

  DetectorMessenger* detectorMessenger;
  RunAction*         run; 

  PrimaryGeneratorAction* gen;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

