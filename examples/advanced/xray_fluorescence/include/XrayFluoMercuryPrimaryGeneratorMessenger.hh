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
// $Id: XrayFluoMercuryPrimaryGeneratorMessenger.hh
// GEANT4 tag $Name:
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

