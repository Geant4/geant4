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
// $Id: Tst52PhysicsListMessenger.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author:Susanna Guatelli (guatelli@ge.infn.it)
//

#ifndef Tst52PhysicsListMessenger_h
#define Tst52PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst52PhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class Tst52PhysicsListMessenger: public G4UImessenger
{
  
public:

  Tst52PhysicsListMessenger(Tst52PhysicsList*);
  ~Tst52PhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  Tst52PhysicsList*          Tst52List;

  G4UIdirectory* EnDir;
   G4UIcmdWithAString*  physicsListCmd;

  G4UIcmdWithADoubleAndUnit* cutECmd; 
  G4UIcmdWithADoubleAndUnit* cutELowECmd; 

  G4UIcmdWithADouble* facRangeCmd;
};

#endif








