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
// $Id: Tst14PhysicsListMessenger.hh,v 1.13 2010-08-29 19:50:17 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Created
//
// -------------------------------------------------------------------

// Class description:
// System test for e/gamma, UI for PhysicsList
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef TST14PHYSICSLISTMESSENGER_HH
#define TST14PHYSICSLISTMESSENGER_HH 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst14PhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class Tst14PhysicsListMessenger: public G4UImessenger {

public:
  
  Tst14PhysicsListMessenger(Tst14PhysicsList* physList);
  
  ~Tst14PhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  Tst14PhysicsList*          physicsList;    // not owned pointer
  
  G4UIdirectory* lowEnDir;

  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;

  G4UIcmdWithADoubleAndUnit* cutGCmd;
  G4UIcmdWithADoubleAndUnit* cutECmd;

  G4UIcmdWithAString*        physicsListCmd;

};

#endif








