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

// ====================================================================
//
//   G4HepMCGenericPythiaMessenger.hh
//   $Id: G4HepMCPythiaMessenger.hh,v 1.1 2002-04-29 20:43:54 asaim Exp $
//
// ====================================================================
#ifndef G4_HEPMC_PYTHIA_MESSENGER_H
#define G4_HEPMC_PYTHIA_MESSENGER_H

#include "G4UImessenger.hh"

class G4HepMCPythiaInterface;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class G4HepMCPythiaMessenger : public G4UImessenger {
private:
  G4HepMCPythiaInterface* gen;

  G4UIdirectory*           dir;
  G4UIcmdWithAnInteger*    verbose;
  G4UIcmdWithAnInteger*    mpylist;
  G4UIcmdWithoutParameter* print;
  G4UIcommand*             cpyinit;
  G4UIcmdWithAnInteger*    cpystat;
  G4UIcommand*             cpygive;
  G4UIcommand*             setUserParameters;
  G4UIcmdWithAnInteger*    setSeed;
  G4UIcommand*             cpyrget;
  G4UIcommand*             cpyrset;
  G4UIcmdWithAString*      printRandomStatus;
  
public:
  G4HepMCPythiaMessenger(G4HepMCPythiaInterface* agen);
  ~G4HepMCPythiaMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
};

#endif
