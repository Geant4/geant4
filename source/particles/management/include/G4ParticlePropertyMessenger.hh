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
// G4ParticlePropertyMessenger
//
// Class description:
//
// This is a messenger class to exchange information between
// G4ParticleDefinition and UI.
//
// /particle/property/   Particle Table control commands.
//   Commands : 
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     Verbose * Set Verbose level

// Author: H.Kurashige, 13 June 1997 - 1st version created
// --------------------------------------------------------------------
#ifndef G4ParticlePropertyMessenger_hh
#define G4ParticlePropertyMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleTable;
class G4ParticleDefinition;
class G4DecayTableMessenger;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger; 

class G4ParticlePropertyMessenger : public G4UImessenger
{
  public:

    G4ParticlePropertyMessenger(G4ParticleTable* pTable = nullptr);
    virtual ~G4ParticlePropertyMessenger();

    G4ParticlePropertyMessenger(const G4ParticlePropertyMessenger&) = delete;
    G4ParticlePropertyMessenger& operator=(const G4ParticlePropertyMessenger&) = delete;

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand* command);

  private:

    G4ParticleTable* theParticleTable = nullptr;

    G4UIdirectory*             thisDirectory = nullptr;
    G4UIcmdWithoutParameter*   dumpCmd = nullptr;
    G4UIcmdWithABool*          stableCmd = nullptr; 
    G4UIcmdWithAnInteger*      verboseCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* lifetimeCmd = nullptr; 
 
    G4DecayTableMessenger* fDecayTableMessenger = nullptr;
};

#endif




