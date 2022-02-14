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
// G4DecayTableMessenger
//
// Class description:
//
// Messenger class to interface and exchange information between the
// decay table/decay channel and UI.  
//
// /particle/property/decay/   Decay Table control commands.
//   Commands : 
//     select * Enter index of decay mode.
//     dump * Dump decay mode information.
//     br * Set branching ratio. [0< BR <1.0]

// Author: H.Kurashige, 13 June 1997
// --------------------------------------------------------------------
#ifndef G4DecayTableMessenger_hh
#define G4DecayTableMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleTable;
class G4VDecayChannel;
class G4ParticleDefinition;
class G4DecayTable;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcmdWithADouble;

class G4DecayTableMessenger : public G4UImessenger
{
  public:

    G4DecayTableMessenger(G4ParticleTable* pTable = nullptr);
    virtual ~G4DecayTableMessenger();

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand* command);

    G4DecayTableMessenger(const G4DecayTableMessenger&) = delete;
    G4DecayTableMessenger& operator= (const G4DecayTableMessenger&) = delete;

  private:

    G4ParticleDefinition* SetCurrentParticle();

    G4ParticleTable* theParticleTable = nullptr;
    G4ParticleDefinition* currentParticle = nullptr;
    G4DecayTable* currentDecayTable = nullptr;
    G4VDecayChannel* currentChannel = nullptr;

    G4UIdirectory* thisDirectory = nullptr;
    G4UIcmdWithoutParameter* dumpCmd = nullptr;
    G4UIcmdWithAnInteger* selectCmd = nullptr;
    G4UIcmdWithADouble* brCmd = nullptr; 

    G4int idxCurrentChannel = -1;
};

#endif
