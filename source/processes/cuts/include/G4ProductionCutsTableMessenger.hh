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
// G4ProductionCutsTableMessenger
//
// Class description:
//
// This is a messenger class to interface to information exchange
// between G4ProductionCutsTable and UI.
// 
// List of Directory and Commands:
//      
//  /run/particle/   Particle control commands
//  Commands : 
//    SetCuts  * Set default cut value
//    dumpList * Dump List of particles in G4VUserPhysicsList.
//    verbose  * Set the Verbose level of G4VUserPhysicsList

// Author: H.Kurashige, 02 March 2008 - First version
// --------------------------------------------------------------------
#ifndef G4ProcductionCutsTableMessenger_hh
#define G4ProcductionCutsTableMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ProductionCutsTable;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString; 
class G4UIcommand;

class G4ProductionCutsTableMessenger : public G4UImessenger
{
  public:

    G4ProductionCutsTableMessenger(G4ProductionCutsTable* pTable);
    virtual ~G4ProductionCutsTableMessenger();

    G4ProductionCutsTableMessenger(const G4ProductionCutsTableMessenger&) = delete;
    G4ProductionCutsTableMessenger& operator=(const G4ProductionCutsTableMessenger&) = delete;
      // Copy contructor and assignment operator not allowed
    
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand* command);

  protected:

    G4ProductionCutsTable* theCutsTable = nullptr;
    
  private:

    G4ProductionCutsTableMessenger() {}
      // Hidden default constructor

    G4UIdirectory*             theDirectory = nullptr;
    G4UIcmdWithAnInteger*      verboseCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* setLowEdgeCmd = nullptr; 
    G4UIcmdWithADoubleAndUnit* setHighEdgeCmd = nullptr; 
    G4UIcmdWithADoubleAndUnit* setMaxEnergyCutCmd = nullptr; 
    G4UIcmdWithoutParameter*   dumpCmd = nullptr;
}; 

#endif
