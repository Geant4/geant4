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
//   H01PythiaMessenger.hh
//   $Id: H01PythiaMessenger.hh,v 1.1 2002-11-19 10:30:55 murakami Exp $
//
// ====================================================================
#ifndef H01_PYTHIA_MESSENGER_H
#define H01_PYTHIA_MESSENGER_H

#include "G4UImessenger.hh"

class H01PythiaInterface;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class H01PythiaMessenger : public G4UImessenger {
private:
  H01PythiaInterface* gen;

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
  H01PythiaMessenger(H01PythiaInterface* agen);
  ~H01PythiaMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
};

#endif
