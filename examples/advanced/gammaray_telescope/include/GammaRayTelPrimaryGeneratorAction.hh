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
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by G.Santin, F.Longo & R.Giannitrapani (30 nov 2000)
//
// ************************************************************



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPrimaryGeneratorAction_h
#define GammaRayTelPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class GammaRayTelDetectorConstruction;
class GammaRayTelPrimaryGeneratorMessenger;
class G4GeneralParticleSource;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  GammaRayTelPrimaryGeneratorAction();    
  ~GammaRayTelPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetSourceType(G4int val) { nSourceType = val;}
  void SetSpectrumType(G4int val) { nSpectrumType = val;}
  void SetVertexRadius(G4double val) { dVertexRadius = val;}
  void SetSourceGen(G4bool val) { sourceGun = val;}
  
private:
  G4ParticleGun*                particleGun;
  G4GeneralParticleSource*      particleSource;	  
  const GammaRayTelDetectorConstruction*    GammaRayTelDetector;  
  GammaRayTelPrimaryGeneratorMessenger* gunMessenger; 
  G4String                      rndmFlag;    //flag for a random impact point
  G4int                         nSourceType;
  G4double                      dVertexRadius;
  G4int                         nSpectrumType;
  G4bool                        sourceGun; // false for GeneralParticleSource
  
};

#endif



