// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XrayFluoPrimaryGeneratorAction.hh,v 1.1 2001-11-27 14:59:32 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoPrimaryGeneratorAction_h
#define XrayFluoPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class XrayFluoDetectorConstruction;
class XrayFluoPrimaryGeneratorMessenger;
class XrayFluoRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    XrayFluoPrimaryGeneratorAction(XrayFluoDetectorConstruction*); 
   
   ~XrayFluoPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetRndmVert (G4String val) { beam = val;}
  void SetSpectrum (G4String val) { spectrum= val  ;}
  void SetIsoVert  (G4String val) { isoVert = val  ;}

private:
  G4ParticleGun*                particleGun;	  //pointer a to G4 service class
  XrayFluoDetectorConstruction*    XrayFluoDetector;  //pointer to the geometry

  XrayFluoPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  XrayFluoRunAction*  runManager;


  G4String                      rndmFlag;   //flag for a random impact point 
  G4String                      beam;
  G4String                      spectrum;
  G4String                      isoVert;
};

#endif


