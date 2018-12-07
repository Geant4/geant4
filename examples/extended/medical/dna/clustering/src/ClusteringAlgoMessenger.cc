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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file ClusteringAlgoMessenger.cc
/// \brief Implementation of the ClusteringAlgoMessenger class

#include "ClusteringAlgoMessenger.hh"

#include "ClusteringAlgo.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgoMessenger::ClusteringAlgoMessenger(ClusteringAlgo* pClustAlgo)
:G4UImessenger(),fpClusteringAlgo(pClustAlgo)
{
  fpAppliDir = new G4UIdirectory("/clustering/");
  fpAppliDir->SetGuidance("commands specific to this example");

  fpMinPtsCmd = new G4UIcmdWithAnInteger("/clustering/algo/setMinPts",this);
  fpMinPtsCmd->SetGuidance("Minimal number of points to create a cluster");
  fpMinPtsCmd->SetParameterName("MinPts",false);
  fpMinPtsCmd->SetRange("MinPts>0");
  fpMinPtsCmd->AvailableForStates(G4State_Idle);

  fpProbCmd = new G4UIcmdWithADouble("/clustering/algo/setSelectionProb",this);
  fpProbCmd->SetGuidance("Probability to select potential "
      "damage according to the geometry");
  fpProbCmd->SetParameterName("Prob",false);
  fpProbCmd->SetRange("Prob>0");
  fpProbCmd->AvailableForStates(G4State_Idle);

  fpEpsCmd = new G4UIcmdWithADoubleAndUnit("/clustering/algo/setEps",this);
  fpEpsCmd->SetGuidance("Maximal distance between points to create a cluster");
  fpEpsCmd->SetParameterName("Eps",false);
  fpEpsCmd->SetRange("Eps>0");
  fpEpsCmd->AvailableForStates(G4State_Idle);

  fpEminCmd = new G4UIcmdWithADoubleAndUnit("/clustering/algo/setEmin",this);
  fpEminCmd->SetGuidance("Energy to have a probability "
      "to create a strand break = 0");
  fpEminCmd->SetParameterName("Emin",false);
  fpEminCmd->SetRange("Emin>=0");
  fpEminCmd->AvailableForStates(G4State_Idle);

  fpEmaxCmd = new G4UIcmdWithADoubleAndUnit("/clustering/algo/setEmax",this);
  fpEmaxCmd->SetGuidance("Energy to have a probability "
      "to create a strand break = 1");
  fpEmaxCmd->SetParameterName("Emax",false);
  fpEmaxCmd->SetRange("Emax>=0");
  fpEmaxCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgoMessenger::~ClusteringAlgoMessenger()
{
  delete fpMinPtsCmd;
  delete fpProbCmd;
  delete fpEpsCmd;
  delete fpEminCmd;
  delete fpEmaxCmd;
  delete fpAppliDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgoMessenger::SetNewValue(G4UIcommand* pCommand,
    G4String pNewValue)
{ 
  if(pCommand == fpMinPtsCmd)
  {
    fpClusteringAlgo->SetMinPts(fpMinPtsCmd->GetNewIntValue(pNewValue));
  }
  if(pCommand == fpProbCmd)
  {
    fpClusteringAlgo->SetSPointsProb(fpProbCmd->GetNewDoubleValue(pNewValue));
  }
  if(pCommand == fpEpsCmd)
  {
    fpClusteringAlgo->SetEps(fpEpsCmd->GetNewDoubleValue(pNewValue));
  }
  if(pCommand == fpEminCmd)
  {
    fpClusteringAlgo->SetEMinDamage(fpEminCmd->GetNewDoubleValue(pNewValue));
  }
  if(pCommand == fpEmaxCmd)
  {
    fpClusteringAlgo->SetEMaxDamage(fpEmaxCmd->GetNewDoubleValue(pNewValue));
  }
}
