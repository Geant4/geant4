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
// G4ProcessManagerMessenger
//
// Class description:
//
// This is a messenger class to exchange information between
// the Process Manager and UI.
//
// /particle/process/   Process Manager control commands
// Commands : 
//     dump * dump process manager information
//     verbose * Set Verbose Level for Process Manager and/or process
//     activate * Activate process  
//     inactivate * Inctivate process  

// Author: H.Kurashige, 13 June 1997
//---------------------------------------------------------------------
#ifndef G4ProcessManagerMessenger_hh
#define G4ProcessManagerMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4ParticleTable;
class G4ProcessManager;
class G4ProcessVector;
class G4VProcess;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcommand;

class G4ProcessManagerMessenger : public G4UImessenger
{
  public:

    G4ProcessManagerMessenger(G4ParticleTable* pTable = nullptr);
      // Constructor

    virtual ~G4ProcessManagerMessenger();
      // Destructor 
 
    G4ProcessManagerMessenger(const G4ProcessManagerMessenger&) = delete;
    G4ProcessManagerMessenger& operator=(const G4ProcessManagerMessenger&) = delete;
      // Copy contructor and assignment operator not allowed

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
      // Set new value for command string

    virtual G4String GetCurrentValue(G4UIcommand* command);
      // Get current value for command string
  
  private:

    const G4ParticleDefinition* SetCurrentParticle();
      // Set particle currently concerned 
    
    G4ParticleTable* theParticleTable = nullptr;
    const G4ParticleDefinition* currentParticle = nullptr;
    G4VProcess* currentProcess = nullptr;
    G4ProcessManager* theManager = nullptr;
    G4ProcessVector* theProcessList = nullptr;

    G4UIdirectory*        thisDirectory = nullptr;
    G4UIcmdWithAnInteger* dumpCmd = nullptr;
    G4UIcommand*          verboseCmd = nullptr;
    G4UIcmdWithAnInteger* activateCmd = nullptr;
    G4UIcmdWithAnInteger* inactivateCmd = nullptr;
};

#endif
