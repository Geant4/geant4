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
// GEANT4 tag $Name: 
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


