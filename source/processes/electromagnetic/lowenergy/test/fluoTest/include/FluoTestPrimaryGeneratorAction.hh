//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestPrimaryGeneratorAction_h
#define FluoTestPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
//class G4ParticleGun;
class G4Event;
class FluoTestDetectorConstruction;
class FluoTestPrimaryGeneratorMessenger;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
   FluoTestPrimaryGeneratorAction(FluoTestDetectorConstruction*);    
   ~FluoTestPrimaryGeneratorAction();

  public:
    
    void GeneratePrimaries(G4Event*);
  // void SetRndmFlag(G4String val) { rndmFlag = val;}
  // void SetRandomizePrimary (G4String val){ randomizePrimary = val;}   
  private:
   G4GeneralParticleSource*                particleGun;	  //pointer a to G4 service class
    FluoTestDetectorConstruction*    Detector;  //pointer to the geometry
    
  // FluoTestPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  //  G4String                      rndmFlag;	  //flag for a random impact point       
  //G4String                      randomizePrimary; //flag for a random energy of the 
                                                   //particle
   
  /*G4double momentum;
    G4double sigmaMomentum;
    G4double sigmaAngle;
 
public:
 
    inline void SetMomentum(G4double val) { momentum = val; }
    inline G4double GetMomentum() const { return momentum; }
    void SetSigmaMomentum(G4double);
    void SetSigmaAngle(G4double);*/  
};

#endif


