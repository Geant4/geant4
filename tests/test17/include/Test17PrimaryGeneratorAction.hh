// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Test17PrimaryGeneratorAction.hh,v 1.1 2000-05-26 06:34:28 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17PrimaryGeneratorAction_h
#define Test17PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Test17DetectorConstruction;
class Test17PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Test17PrimaryGeneratorAction(Test17DetectorConstruction*);    
   ~Test17PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();   
    void GeneratePrimaries(G4Event*);
    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;
    Test17DetectorConstruction*      Test17Detector;
    
    static G4String thePrimaryParticleName;
    
    Test17PrimaryGeneratorMessenger* gunMessenger;     
};

#endif


