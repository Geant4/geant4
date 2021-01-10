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
/// \file eventgenerator/HepMC/HepMCEx02/include/H02PrimaryGeneratorMessenger.hh
/// \brief Definition of the H02PrimaryGeneratorMessenger class
//
//
#ifndef H02_PRIMARY_GENERATOR_MESSENGER_H
#define H02_PRIMARY_GENERATOR_MESSENGER_H

#include "globals.hh"
#include "G4UImessenger.hh"

class H02PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class H02PrimaryGeneratorMessenger : public G4UImessenger {
public:
  H02PrimaryGeneratorMessenger(H02PrimaryGeneratorAction* genaction);
  ~H02PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand* command, G4String newValues);
  virtual G4String GetCurrentValue(G4UIcommand* command);

private:
  H02PrimaryGeneratorAction* fPrimaryAction;

  G4UIdirectory* fDir;
  G4UIcmdWithAString* fSelect;

};

#endif
