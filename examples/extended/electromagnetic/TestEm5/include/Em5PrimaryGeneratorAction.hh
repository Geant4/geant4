// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5PrimaryGeneratorAction.hh,v 1.2 1999-12-15 14:49:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5PrimaryGeneratorAction_h
#define Em5PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Em5DetectorConstruction;
class Em5PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Em5PrimaryGeneratorAction(Em5DetectorConstruction*);    
   ~Em5PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();   
    void GeneratePrimaries(G4Event*);
    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;
    Em5DetectorConstruction*      Em5Detector;
    
    static G4String thePrimaryParticleName;
    
    Em5PrimaryGeneratorMessenger* gunMessenger;     
};

#endif


