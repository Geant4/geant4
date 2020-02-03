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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file EventActionMessenger.cc
/// \brief Implementation of the EventActionMessenger class

#include "EventActionMessenger.hh"

#include "EventAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventActionMessenger::EventActionMessenger(EventAction* EvAct)
:G4UImessenger(),fpEventAction(EvAct)
{
  fpPDBDir = new G4UIdirectory("/PDB4DNA/");
  fpPDBDir->SetGuidance("commands specific to this example");

  fpThresEdepCmd = new G4UIcmdWithADoubleAndUnit(
      "/PDB4DNA/event/setEnergyThres",
      this);
  fpThresEdepCmd->SetGuidance("Set energy threshold for SSB");
  fpThresEdepCmd->SetParameterName("EnergyThres",false);
  fpThresEdepCmd->SetRange("EnergyThres>0");
  fpThresEdepCmd->AvailableForStates(G4State_Idle);

  fpThresDistCmd = new G4UIcmdWithAnInteger("/PDB4DNA/event/setDistanceThres",
                                            this);
  fpThresDistCmd->SetGuidance("Set distance threshold for DSB");
  fpThresDistCmd->SetParameterName("DistanceThres",false);
  fpThresDistCmd->SetRange("DistanceThres>0");
  fpThresDistCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventActionMessenger::~EventActionMessenger()
{
  delete fpThresEdepCmd;
  delete fpThresDistCmd;
  delete fpPDBDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command==fpThresEdepCmd)
  {
    fpEventAction->SetEnergyThresForSSB(fpThresEdepCmd->GetNewDoubleValue(
        newValue));
  }

  if(command==fpThresDistCmd)
  {
    fpEventAction->SetDistanceThresForDSB(fpThresDistCmd->GetNewIntValue(
        newValue));
  }

}

