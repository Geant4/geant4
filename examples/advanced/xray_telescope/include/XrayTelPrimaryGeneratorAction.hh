// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorAction.hh                *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

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

  G4ThreeVector GetRandomRingPosition (G4double minRad, G4double maxRad);   
  // returns a random position vector in a ring defined by minRad, maxRad
  G4ThreeVector GetRandomDirection(G4double maxT);
  // returns a random direction vector in the space
  HepJamesRandom* engine;
};

#endif


