// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationMessenger.hh,v 1.2 1999-04-14 14:25:24 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4UIdirectory      *fFSDirectory;
  G4UIcmdWithAString *fListEnvelopesCmd;
  G4UIcmdWithAString *fListModelsCmd;
  G4UIcmdWithAString *fListIsApplicableCmd;
  G4UIcmdWithAString *fActivateModel;
  G4UIcmdWithAString *fInActivateModel;
};

#endif


