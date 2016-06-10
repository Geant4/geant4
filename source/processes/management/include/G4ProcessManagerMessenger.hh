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
//
// $Id: G4ProcessManagerMessenger.hh 71231 2013-06-12 13:06:28Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ProcessManagerMessenger.hh
//
// Class Description:
//    This is a messenger class to interface to exchange information
//    between ProcessManager and UI.
//-
//  /particle/process/   Process Manager control commands.
//   Commands : 
//     dump * dump process manager information.
//     Verbose * Set Verbose Level for Process Manager and/or process
//     Activate * Activate process  
//     inactivate * Inctivate process  
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov. 1997  H. Kurashige   : fixed bugs 
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//
//---------------------------------------------------------------

#ifndef G4ProcessManagerMessenger_h
#define G4ProcessManagerMessenger_h 1

class G4ParticleDefinition;
class G4ParticleTable;
class G4ProcessManager;
class G4ProcessVector;
class G4VProcess;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ProcessManagerMessenger: public G4UImessenger
{
  public:
    G4ProcessManagerMessenger(G4ParticleTable* pTable = 0);
    // constructor

    virtual ~G4ProcessManagerMessenger();
    // destructor 
 
public: // with description
     virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    // set new value for command string

    virtual G4String GetCurrentValue(G4UIcommand * command);
    // get current value for command string
  
  private:
    G4ProcessManagerMessenger(const G4ProcessManagerMessenger&)
      : G4UImessenger() {};
    // hide copy constructor as private

  private:
    G4ParticleDefinition* SetCurrentParticle();
    // set particle currently concerned 
    
  private:
    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;
    G4VProcess* currentProcess;
    G4ProcessManager*     theManager;
    G4ProcessVector*      theProcessList;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithAnInteger *      dumpCmd;
    G4UIcommand          *      verboseCmd;
    G4UIcmdWithAnInteger *      activateCmd;
    G4UIcmdWithAnInteger *      inactivateCmd;
};

#endif

