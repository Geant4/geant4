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
// $Id: RemSimPhysicsListMessenger.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef RemSimPhysicsListMessenger_h
#define RemSimPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RemSimPhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class RemSimPhysicsListMessenger: public G4UImessenger
{
  
public:

  RemSimPhysicsListMessenger(RemSimPhysicsList*);
  ~RemSimPhysicsListMessenger();
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  RemSimPhysicsList*          RemSimList;
  G4UIdirectory* EnDir;
   G4UIcmdWithAString*  physicsListCmd;
  G4UIcmdWithADoubleAndUnit* cutECmd;
};

#endif








