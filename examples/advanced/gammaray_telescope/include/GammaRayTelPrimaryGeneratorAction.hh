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

  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4 service class
    GammaRayTelDetectorConstruction*    GammaRayTelDetector;  //pointer to the geometry
    
    GammaRayTelPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;	  //flag for a random impact point       
};

#endif


