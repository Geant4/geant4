// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorAction.hh,v 1.7 2001-10-15 09:21:07 guardi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestPrimaryGeneratorAction_h
#define FluoTestPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class FluoTestDetectorConstruction;
class FluoTestPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    FluoTestPrimaryGeneratorAction(FluoTestDetectorConstruction*);    
   ~FluoTestPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}

  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4 service class
    FluoTestDetectorConstruction*    FluoTestDetector;  //pointer to the geometry
    
    FluoTestPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;	  //flag for a random impact point       
};

#endif


