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
// $Id: G4VisCommandsMultithreading.hh 66264 2012-12-14 10:17:44Z allison $

// /vis/multithreading commands - John Allison  29th September 2015

#ifdef G4MULTITHREADED

#ifndef G4VISCOMMANDSMULTITHREADING_HH
#define G4VISCOMMANDSMULTITHREADING_HH

#include "G4VVisCommand.hh"

class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class G4VisCommandMultithreadingActionOnEventQueueFull: public G4VVisCommand {
public:
  G4VisCommandMultithreadingActionOnEventQueueFull();
  virtual ~G4VisCommandMultithreadingActionOnEventQueueFull();
  G4String GetCurrentValue(G4UIcommand* command);
  void SetNewValue(G4UIcommand* command, G4String newValue);
private:
  G4VisCommandMultithreadingActionOnEventQueueFull
  (const G4VisCommandMultithreadingActionOnEventQueueFull&);
  G4VisCommandMultithreadingActionOnEventQueueFull& operator=
  (const G4VisCommandMultithreadingActionOnEventQueueFull&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandMultithreadingMaxEventQueueSize: public G4VVisCommand {
public:
  G4VisCommandMultithreadingMaxEventQueueSize();
  virtual ~G4VisCommandMultithreadingMaxEventQueueSize();
  G4String GetCurrentValue(G4UIcommand* command);
  void SetNewValue(G4UIcommand* command, G4String newValue);
private:
  G4VisCommandMultithreadingMaxEventQueueSize
  (const G4VisCommandMultithreadingMaxEventQueueSize&);
  G4VisCommandMultithreadingMaxEventQueueSize& operator=
  (const G4VisCommandMultithreadingMaxEventQueueSize&);
  G4UIcmdWithAnInteger* fpCommand;
};

#endif

#endif
