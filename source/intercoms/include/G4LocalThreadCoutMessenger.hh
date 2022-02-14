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
// G4LocalThreadCoutMessenger
//
// Class description
//
// This class is the messenger for handling cout/cerr of a local thread

// Author: M.Asai, 2013
// --------------------------------------------------------------------
#ifndef G4LocalThreadCoutMessenger_hh
#define G4LocalThreadCoutMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class G4LocalThreadCoutMessenger : public G4UImessenger
{
  public:

    G4LocalThreadCoutMessenger();
    ~G4LocalThreadCoutMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:

    G4UIdirectory* coutDir = nullptr;
    G4UIcommand* coutFileNameCmd = nullptr;
    G4UIcommand* cerrFileNameCmd = nullptr;
    G4UIcmdWithABool* bufferCoutCmd = nullptr;
    G4UIcmdWithAString* prefixCmd = nullptr;
    G4UIcmdWithAnInteger* ignoreCmd = nullptr;
    G4UIcmdWithABool* ignoreInitCmd = nullptr;
};

#endif
