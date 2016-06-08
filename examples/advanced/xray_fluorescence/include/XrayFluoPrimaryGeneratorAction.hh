//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoPrimaryGeneratorAction.hh
// GEANT4 tag $Name:  xray_fluo-V04-01-03
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

  //method to set a random impact point on the sample
  void SetRndmFlag(G4String val) { rndmFlag = val;}

  //method to choose a circular non-point source
  void SetRndmVert (G4String val) { beam = val;}

  //set the flag for shooting particles according to certain spectra 
  void SetSpectrum (G4String val) { spectrum= val  ;}

 //set the flag for shooting particles from an isotropic source
  void SetIsoVert  (G4String val) { isoVert = val  ;}

private:
//pointer a to G4 service class
  G4ParticleGun*                particleGun;	  

 //pointer to the geometry
  XrayFluoDetectorConstruction*    XrayFluoDetector;  

  //messenger of this class
  XrayFluoPrimaryGeneratorMessenger* gunMessenger; 
  XrayFluoRunAction*  runManager;

 //flag for a random impact point 
  G4String                      rndmFlag;   

  //flag for a circular non-point source
  G4String                      beam;

 //flag  for shooting particles according to certain spectra 
  G4String                      spectrum;

 //flag  for shooting particles from an isotropic source
  G4String                      isoVert;
};

#endif


