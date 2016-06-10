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
// $Id: G4FastSimulationMessenger.hh 68056 2013-03-13 14:44:48Z gcosmo $
//
// 
//	GEANT 4 class header file 
//
//      This is a messenger class for G4FastSimulation.
//      Implemented commands are following;
//
//  Commands : 
//    BeamOn *              Start a Run.
// 
//	History
//        first version  by P.Mora de Freitas & M.Verderi 
// ------------------------------------------------------------

#ifndef G4FastSimulationMessenger_h
#define G4FastSimulationMessenger_h 1

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4GlobalFastSimulationManager.hh"

class G4FastSimulationMessenger: public G4UImessenger
{
public:
  G4FastSimulationMessenger(G4GlobalFastSimulationManager* theGFSM);
  virtual ~G4FastSimulationMessenger();
  
public:
  void SetNewValue(G4UIcommand * command,G4String newValues);
  
private:
  G4GlobalFastSimulationManager* fGlobalFastSimulationManager;
  
  //commands
  G4UIdirectory*                   fFSDirectory;
  G4UIcmdWithoutParameter*        fShowSetupCmd;
  G4UIcmdWithAString*         fListEnvelopesCmd;
  G4UIcmdWithAString*            fListModelsCmd;
  G4UIcmdWithAString*      fListIsApplicableCmd;
  G4UIcmdWithAString*            fActivateModel;
  G4UIcmdWithAString*          fInActivateModel;
};

#endif


