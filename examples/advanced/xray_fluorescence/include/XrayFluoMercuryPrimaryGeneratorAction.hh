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
// $Id: XrayFluoPlanePrimaryGeneratorAction.hh
// GEANT4 tag $Name:  
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef XrayFluoMercuryPrimaryGeneratorAction_h
#define XrayFluoMercuryPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class XrayFluoMercuryDetectorConstruction;
class XrayFluoMercuryPrimaryGeneratorMessenger;
class XrayFluoRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoMercuryPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    XrayFluoMercuryPrimaryGeneratorAction(XrayFluoMercuryDetectorConstruction*); 
   
   ~XrayFluoMercuryPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event*);

  // method to set the global illumination of Mercury surface 
  void SetGlobalFlag (G4bool val) {globalFlag = val;};


  //method to set a random impact point on the sample
  //void SetRndmFlag(G4String val) { rndmFlag = val;}

  //method to choose a circular non-point-like source
  //void SetRndmVert (G4String val) { beam = val;};

  //set the flag for shooting particles according to certain spectra 
  void SetSpectrum (G4String val) { spectrum = val;};

  //set the flag for shooting particles from an isotropic source
  //void SetIsoVert  (G4String val) { isoVert = val  ;};

private:
//pointer a to G4 service class
  G4ParticleGun*                particleGun;	  

  //pointer to the geometry
  XrayFluoMercuryDetectorConstruction*    XrayFluoDetector;  
  
  //messenger of this class
  XrayFluoMercuryPrimaryGeneratorMessenger* gunMessenger; 
  XrayFluoRunAction*  runManager;
  
  //flag for setting the global illumination of Mercury surface

  G4bool                        globalFlag;

  //flag for considering the Sun as a Point-like source
  //not to be used if spacecraft latitude is near zero
  //because the spacecraft itself shadowes the planet surface.

  //G4bool                        pointLikeFlag;

  //flag for a random impact point 
  //  G4String                      rndmFlag;   
  
  //flag for a circular non-point source
  //  G4String                      beam;

 //flag  for shooting particles according to certain spectra 
  G4String                      spectrum;

 //flag  for shooting particles from an isotropic source
  //  G4String                      isoVert;
};

#endif
