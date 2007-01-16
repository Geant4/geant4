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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
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
  G4UIcmdWithoutParameter*       cleanCmd;

  G4String                       bodypart;
  G4bool                         bps;

};

#endif

