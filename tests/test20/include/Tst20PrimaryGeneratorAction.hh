// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20PrimaryGeneratorAction.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20PrimaryGeneratorAction  ------
// ************************************************************



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20PrimaryGeneratorAction_h
#define Tst20PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Tst20DetectorConstruction;
class Tst20PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst20PrimaryGeneratorAction(Tst20DetectorConstruction*);    
  ~Tst20PrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetSourceType(G4int val) { nSourceType = val;}
  void SetSpectrumType(G4int val) { nSpectrumType = val;}
  void SetVertexRadius(G4double val) { dVertexRadius = val;}
  
private:
  G4ParticleGun*                particleGun;	  
  Tst20DetectorConstruction*    Tst20Detector;  
  Tst20PrimaryGeneratorMessenger* gunMessenger; 
  G4String                      rndmFlag;    //flag for a random impact point
  G4int                         nSourceType;
  G4double                      dVertexRadius;
  G4int                         nSpectrumType;
};

#endif



