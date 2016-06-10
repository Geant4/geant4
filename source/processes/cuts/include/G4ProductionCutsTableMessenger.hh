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
// $Id: G4ProductionCutsTableMessenger.hh 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
//---------------------------------------------------------------
//
//  G4ProcductionCutsTableMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between ProductionCutsTable and UI.
// --
//  the List of Directory and Commands
// -       
//  /run/particle/   Paricle control commands.
//   Commands : 
//    SetCuts *  Set default cut value
//    dumpList * Dump List of particles in G4VUserPhysicsList.
//    verbose * Set the Verbose level of G4VUserPhysicsList.
// ------------------------------------------------------------
//	History
//        first version                02 Mar. 2008 by H.Kurashige 
// ------------------------------------------------------------

#ifndef G4ProcductionCutsTableMessenger_h
#define G4ProcductionCutsTableMessenger_h 1

class G4ProductionCutsTable;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString; 
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ProductionCutsTableMessenger: public G4UImessenger
{
  private:
  // hide default constructor
    G4ProductionCutsTableMessenger(){}

  public:
    G4ProductionCutsTableMessenger(G4ProductionCutsTable* pTable);
    virtual ~G4ProductionCutsTableMessenger();
    
public: // with description
    virtual  void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual  G4String GetCurrentValue(G4UIcommand * command);

  protected:
    G4ProductionCutsTable* theCutsTable;
    
  private: //commands
    G4UIdirectory *             theDirectory;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithADoubleAndUnit * setLowEdgeCmd; 
    G4UIcmdWithADoubleAndUnit * setHighEdgeCmd; 
    G4UIcmdWithADoubleAndUnit * setMaxEnergyCutCmd; 
    G4UIcmdWithoutParameter *   dumpCmd;
 
}; 

#endif


