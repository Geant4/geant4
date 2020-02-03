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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  19 Sep 03  Alfonso Mantero created  
//
// -------------------------------------------------------------------


#ifndef XrayFluoMercuryPrimaryGeneratorMessenger_h
#define XrayFluoMercuryPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
//#include "XrayFluoPlanePrimaryGeneratorAction.hh"

class XrayFluoMercuryPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoMercuryPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  XrayFluoMercuryPrimaryGeneratorMessenger(XrayFluoMercuryPrimaryGeneratorAction*);
  ~XrayFluoMercuryPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  XrayFluoMercuryPrimaryGeneratorAction* XrayFluoAction; 

//   //command to set a random impact point
//   G4UIcmdWithAString*          RndmCmd;

//  //command to choose a plane circular source
//   G4UIcmdWithAString*          RndmVert;

  //command to shoot particles according to certain spectra
  G4UIcmdWithAString* spectrum;

  //command to shoot particles to all Mercury illuminated surface
  G4UIcmdWithABool* globalCmd;

};

#endif

