//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef XrayFluoPrimaryGeneratorAction_h
#define XrayFluoPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4ParticleGun;
class G4Event;
class XrayFluoDetectorConstruction;
class XrayFluoPrimaryGeneratorMessenger;
class XrayFluoRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    XrayFluoPrimaryGeneratorAction(const XrayFluoDetectorConstruction*); 
   
   ~XrayFluoPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event*);

  //method to set a random impact point on the sample
  void SetRndmFlag(G4String val) { rndmFlag = val;}

  //method to choose a circular non-point source
  void SetRndmVert (G4String val) { beam = val;}

  //set the flag for shooting particles according to certain spectra 
  void SetSpectrum (G4String val) { spectrum= val  ;}

 //set the flag for shooting particles from an isotropic source
  void SetIsoVert  (G4String val) { isoVert = val  ;}

  //load data created in previous runs 
  void ActivatePhaseSpace(G4String);

  //load Raileigh scattered photon from previously loaded data
  void SetRayleighFlag (G4bool);


private:
//pointer a to G4 service class
  G4ParticleGun*                particleGun;	  

  //pointer to the geometry
  const XrayFluoDetectorConstruction*    XrayFluoDetector;  
  
  //pointer to the run manager (master)
  const XrayFluoRunAction* runAction;

  //messenger of this class
  XrayFluoPrimaryGeneratorMessenger* gunMessenger; 
  
  //flag for a random impact point 
  G4String                      rndmFlag;   
  
  //flag for a circular non-point source
  G4String                      beam;

 //flag  for shooting particles according to certain spectra 
  G4String                      spectrum;

 //flag  for shooting particles from an isotropic source
  G4String                      isoVert;

  //flag for the use of phase-space created with a previous run
  G4bool phaseSpaceGunFlag;

  // the flag to load  particle coming from Rayleigh scattering
  G4bool rayleighFlag; 

  G4ThreeVector detectorPosition;



};

#endif


