//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoPrimaryGeneratorAction_h
#define XrayFluoPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class XrayFluoDetectorConstruction;
class XrayFluoPrimaryGeneratorMessenger;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  XrayFluoPrimaryGeneratorAction(XrayFluoDetectorConstruction*);    
  ~XrayFluoPrimaryGeneratorAction();
  
public:
  
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetRandomizePrimary (G4String val){ randomizePrimary = val;}   
private:
  G4ParticleGun*                particleGun;	  //pointer a to G4 service class
  XrayFluoDetectorConstruction*    Detector;  //pointer to the geometry
  
  XrayFluoPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  G4String                      rndmFlag;	  //flag for a random impact point       
  G4String                      randomizePrimary; //flag for a random energy of the 
  //particle
  
  G4double momentum;
  G4double sigmaMomentum;
  G4double sigmaAngle;
  
public:
  
  inline void SetMomentum(G4double val) { momentum = val; }
  inline G4double GetMomentum() const { return momentum; }
  void SetSigmaMomentum(G4double);
  void SetSigmaAngle(G4double);
};

#endif


