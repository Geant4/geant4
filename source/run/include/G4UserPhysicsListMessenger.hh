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
// G4UserPhysicsListMessenger
//
// Class description:
//
// This is a messenger class to  allow exchange of information
// between ParticleUserList and UI.
// 
// Directory and list of commands:
// 
// /run/particle/   Particle control commands.
//  Commands :
//   SetCuts  * Set default cut value
//   dumpList * Dump List of particles in G4VUserPhysicsList.
//   verbose  * Set the Verbose level of G4VUserPhysicsList.
//   addProcessManager    * add process manager
//   buildPhysicsTable    * build physics table
//   storePhysicsTable    * store physics table into files
//   retreivePhysicsTable * retrieve physics table from files
//   setStoredInAscii * Switch on/off ascii mode in store/retrieve Physics Table

// Original author: H.Kurashige, 9 January 1998
// --------------------------------------------------------------------
#ifndef G4UserPhysicsListMessenger_hh
#define G4UserPhysicsListMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4VUserPhysicsList;
class G4VUserPhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;

class G4UserPhysicsListMessenger : public G4UImessenger
{
  public:

    G4UserPhysicsListMessenger(G4VUserPhysicsList* pParticleList);
    virtual ~G4UserPhysicsListMessenger();

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand* command);

  protected:

    G4VUserPhysicsList* thePhysicsList = nullptr;

  private:

    G4UserPhysicsListMessenger() {}
      // Hidden default constructor.

    G4UIdirectory* theDirectory = nullptr;
    G4UIcmdWithADoubleAndUnit* setCutCmd = nullptr;
    G4UIcommand* setCutRCmd = nullptr;
    G4UIcommand* setCutForAGivenParticleCmd = nullptr;
    G4UIcmdWithAString* getCutForAGivenParticleCmd = nullptr;
    G4UIcmdWithAnInteger* verboseCmd = nullptr;
    G4UIcmdWithoutParameter* dumpListCmd = nullptr;
    G4UIcmdWithAString* addProcManCmd = nullptr;
    G4UIcmdWithAString* buildPTCmd = nullptr;
    G4UIcmdWithAString* storeCmd = nullptr;
    G4UIcmdWithAString* retrieveCmd = nullptr;
    G4UIcmdWithAnInteger* asciiCmd = nullptr;
    G4UIcommand* applyCutsCmd = nullptr;
    G4UIcmdWithAString* dumpCutValuesCmd = nullptr;
    G4UIcmdWithAnInteger* dumpOrdParamCmd = nullptr;
};

#endif
