//  XrayTelPrimaryGeneratorAction.hh
// 

#ifndef XrayTelPrimaryGeneratorAction_h
#define XrayTelPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Random/RandGeneral.h"

class G4ParticleGun;
class G4Event;
class XrayTelDetectorConstruction;
class XrayTelPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  XrayTelPrimaryGeneratorAction( XrayTelDetectorConstruction* );    
  ~XrayTelPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetErndmFlag(G4String val) { erndmFlag = val;}
  void SetRmin ( G4double r ) { Rmin = r; }
  void SetRmax ( G4double r ) { Rmax = r; }
  void SetTmax ( G4double t ) { Tmax = t; }

private:
  G4ParticleGun* particleGun;	          // pointer a to G4 service class   
  XrayTelDetectorConstruction* XrayTelDetector;  //pointer to the geometry
  XrayTelPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  G4String rndmFlag;	              //flag for a random impact point    
  G4String erndmFlag;	              //flag for a random energy selection    
  G4double Rmin;
  G4double Rmax;
  G4double Tmax;

  G4ThreeVector GetRandomShellPosition (G4double Rmin, G4double Rmax);   
  // returns a random position vector in a spherical shell defined by Rmin, Rmax
  G4ThreeVector GetRandomRingPosition (G4double Rmin, G4double Rmax);   
  // returns a random position vector in a ring defined by Rmin, Rmax
  G4ThreeVector GetRandomDirection();
  G4ThreeVector GetRandomDirection(G4double Tmax);
  // returns a random direction vector in the space
  G4double GetRandomEnergy(G4int ParticleCode);
  // returns the particle energy ( ParticleCode == 1 electrons, ParticleCode == 2 protons )
  // according to the fluxes data spectra
  G4ThreeVector Get2DRandomDirection();
  G4ThreeVector GetRandomPositionOnaCircle ( G4double R );
  HepJamesRandom* engine;
  RandGeneral* ElectronRandomEnergy;
  RandGeneral* ProtonRandomEnergy;
};

#endif


