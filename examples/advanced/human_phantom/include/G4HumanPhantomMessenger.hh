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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//
#ifndef G4HumanPhantomMessenger_h
#define G4HumanPhantomMessenger_h 1

class G4HumanPhantomConstruction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class G4HumanPhantomMessenger: public G4UImessenger
{
public:
  G4HumanPhantomMessenger(G4HumanPhantomConstruction* myUsrPhtm);
  ~G4HumanPhantomMessenger();
    
  void SetNewValue(G4UIcommand* command, G4String newValue);

  void AddBodyPart(G4String);	      // Set Body Parts Sensitivity

private:
  G4HumanPhantomConstruction*           myUserPhantom;

  G4UIdirectory*                 phantomDir;
  G4UIdirectory*                 bpDir;

  G4UIcmdWithAString*            modelCmd; 
  G4UIcmdWithAString*            sexCmd;  
  G4UIcmdWithAString*            bodypartCmd;
  G4UIcmdWithoutParameter*       endCmd;

  G4String                       bodypart;
  G4bool                         bps;

};

#endif

