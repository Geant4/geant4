// Rich advanced example for Geant4
// RichTbPrimaryGeneratorAction.hh for Rich of LHCb
// History:
// Created: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////

#ifndef RichTbPrimaryGeneratorAction_h
#define RichTbPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class RichTbDetectorConstruction;
class RichTbPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RichTbPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    RichTbPrimaryGeneratorAction(RichTbDetectorConstruction*);    
   ~RichTbPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}
    void Setxvertex(G4double x) ;
    void Setyvertex(G4double y) ;
    void Setzvertex(G4double z) ;

    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;	//pointer a to G4 service class
    RichTbDetectorConstruction*      RichTbDetector; //pointer to the geometry
    
    RichTbPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;	//flag for a random impact point       

    static G4String thePrimaryParticleName ;
    G4double xvertex,yvertex,zvertex;
    G4bool vertexdefined;

};

#endif


