// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPrimaryGeneratorAction.hh,v 1.3 2000-12-06 16:53:13 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GammaRayTelPrimaryGeneratorAction(GammaRayTelDetectorConstruction*);    
  ~GammaRayTelPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetSourceType(G4int val) { nSourceType = val;}
  void SetSpectrumType(G4int val) { nSpectrumType = val;}
  void SetVertexRadius(G4double val) { dVertexRadius = val;}
  
private:
  G4ParticleGun*                particleGun;	  
  GammaRayTelDetectorConstruction*    GammaRayTelDetector;  
  GammaRayTelPrimaryGeneratorMessenger* gunMessenger; 
  G4String                      rndmFlag;    //flag for a random impact point
  G4int                         nSourceType;
  G4double                      dVertexRadius;
  G4int                         nSpectrumType;
};

#endif



