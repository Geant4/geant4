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
// $Id: Tst50PhysicsListMessenger.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef Tst50PhysicsListMessenger_h
#define Tst50PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst50PhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst50PhysicsListMessenger: public G4UImessenger
{
  
public:

  Tst50PhysicsListMessenger(Tst50PhysicsList*);
  ~Tst50PhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  Tst50PhysicsList*          Tst50List;

  G4UIdirectory* lowEnDir;
   G4UIcmdWithAString* RangeDir;
 
  G4UIcmdWithADoubleAndUnit* cutECmd;
};

#endif








