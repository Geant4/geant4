//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4FastSimulationMessenger.hh,v 1.4 2001-07-11 10:08:23 gunter Exp $
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


