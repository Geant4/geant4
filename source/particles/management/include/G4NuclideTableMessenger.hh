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
// G4NuclideTableMessenger
//
// Class Description:
//
// This is a messenger class to exchange information between
// G4ParticleDefinition and UI.
//
//   /particle/manage/nuclide   Nuclide Table control commands
//   Commands : 
//     lifetime * Set threshold of half-life.

// Author: T.Koi, SLAC - 11 November 2015
// --------------------------------------------------------------------
#ifndef G4NuclideTableMessenger_hh
#define G4NuclideTableMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4NuclideTable;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

class G4NuclideTableMessenger : public G4UImessenger
{
  public:

    G4NuclideTableMessenger(G4NuclideTable* nuclideTable = nullptr);
    virtual ~G4NuclideTableMessenger();

    G4NuclideTableMessenger(const G4NuclideTableMessenger&) = delete;
    G4NuclideTableMessenger& operator=(const G4NuclideTableMessenger&) = delete;

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

  private:

    G4NuclideTable* theNuclideTable = nullptr;

    G4UIdirectory* thisDirectory = nullptr;
    G4UIcmdWithADoubleAndUnit* halflifeCmd = nullptr; 
    G4UIcmdWithADoubleAndUnit* meanlifeCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* lToleranceCmd = nullptr; 
};

#endif
