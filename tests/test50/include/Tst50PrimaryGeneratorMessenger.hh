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
// $Id: Tst50PrimaryGeneratorMessenger.hh
// GEANT4 tag $Name:  xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef Tst50PrimaryGeneratorMessenger_h
#define Tst50PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Tst50PrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst50PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  Tst50PrimaryGeneratorMessenger(Tst50PrimaryGeneratorAction*);
  ~Tst50PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  Tst50PrimaryGeneratorAction* Tst50Action; 

  //command to set a random impact point
  G4UIcmdWithAString*          RndmCmd;

 //command to choose a plane circular source
  G4UIcmdWithAString*          RndmVert;
 
  //command to shot particles according to certain spectra
  G4UIcmdWithAString*        spectrum;

 //command to shot particles from an isotropic source
  G4UIcmdWithAString*        isoVert;
};

#endif

