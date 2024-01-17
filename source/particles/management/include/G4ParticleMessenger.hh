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
// G4ParticleMessenger
//
// Class description:
//
// This is a messenger class to exchange information between particle
// related classes and UI.
//
// List of Directory and Commands:
//
// G4ParticleMessenger
//  /particle/   Particle control commands.
//   Commands :
//    select * Select particle
//    list * List name of particles.
//    find * find particle by PDG encoding.
//    verbose * Set Verbose level of Particle Table
//
// G4ParticlePropertyMessenger
//  /particle/property/   Particle Table control commands.
//   Commands :
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     verbose * Set Verbose level
//
// G4DecayTableMessenger
//  /particle/property/decay/   Decay Table control commands.
//   Commands :
//     select * Enter index of decay mode.
//     dump * Dump decay mode information.
//     br * Set branching ratio. [0< BR <1.0]

// Author: H. Kurashige, 13 June 1997 - 1st version created
// --------------------------------------------------------------------
#ifndef G4ParticleMessenger_hh
#define G4ParticleMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4ParticleTable;
class G4ParticlePropertyMessenger;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class G4ParticleMessenger : public G4UImessenger
{
  public:
    G4ParticleMessenger(G4ParticleTable* pTable = nullptr);
    ~G4ParticleMessenger() override;

    G4ParticleMessenger(const G4ParticleMessenger&) = delete;
    G4ParticleMessenger& operator=(const G4ParticleMessenger&) = delete;

    void SetNewValue(G4UIcommand* command, G4String newValues) override;
    G4String GetCurrentValue(G4UIcommand* command) override;

  private:
    G4UIdirectory* thisDirectory = nullptr;
    G4UIcmdWithAString* listCmd = nullptr;
    G4UIcmdWithAString* selectCmd = nullptr;
    G4UIcmdWithAnInteger* findCmd = nullptr;
    G4UIcmdWithoutParameter* createAllIonCmd = nullptr;
    G4UIcmdWithoutParameter* createAllIsomerCmd = nullptr;
    G4UIcmdWithAnInteger* verboseCmd = nullptr;

    G4ParticleTable* theParticleTable = nullptr;
    G4ParticleDefinition* currentParticle = nullptr;

    G4ParticlePropertyMessenger* fParticlePropertyMessenger = nullptr;
};

#endif
