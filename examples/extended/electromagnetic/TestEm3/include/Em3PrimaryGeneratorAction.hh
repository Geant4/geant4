// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PrimaryGeneratorAction.hh,v 1.1 1999-10-11 16:55:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3PrimaryGeneratorAction_h
#define Em3PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Em3DetectorConstruction;
class Em3PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Em3PrimaryGeneratorAction(Em3DetectorConstruction*);    
   ~Em3PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic(); 
    void GeneratePrimaries(G4Event*);
    
  private:
    G4ParticleGun*              particleGun;
    Em3DetectorConstruction*    Em3Detector; 
    
    Em3PrimaryGeneratorMessenger* gunMessenger; 
};

#endif


