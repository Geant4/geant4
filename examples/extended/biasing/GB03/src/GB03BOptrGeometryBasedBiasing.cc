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
/// \file GB03BOptrGeometryBasedBiasing.cc
/// \brief Implementation of the GB03BOptrGeometryBasedBiasing class

#include "GB03BOptrGeometryBasedBiasing.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03BOptrGeometryBasedBiasing::GB03BOptrGeometryBasedBiasing()
: G4VBiasingOperator("GB03BOptrGeometryBasedBiasing"),
  fSplittingFactor(2),
  fApplyProbability(1.0)
{
  fSplitAndKillOperation = new GB03BOptnSplitOrKillOnBoundary("splitAndkill");
  
  // -- Define messengers:
  fSplittingFactorMessenger = 
    new G4GenericMessenger(this, "/GB03/biasing/","Biasing control" );
  
  G4GenericMessenger::Command& splittingFactorCmd = 
    fSplittingFactorMessenger->DeclareProperty("setSplittingFactor", fSplittingFactor,
                                "Define the splitting factor." );
  splittingFactorCmd.SetStates(G4State_Idle);
  
  fApplyProbabilityMessenger = 
    new G4GenericMessenger(this, "/GB03/biasing/","Biasing control" );
  
  G4GenericMessenger::Command& applyProbCmd = 
    fApplyProbabilityMessenger->DeclareProperty("setApplyProbability", fApplyProbability,
                                "Define the probability to apply the splitting/killing." );
  applyProbCmd.SetStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03BOptrGeometryBasedBiasing::~GB03BOptrGeometryBasedBiasing()
{
  delete fSplitAndKillOperation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03BOptrGeometryBasedBiasing::StartRun()
{
  fSplitAndKillOperation->SetSplittingFactor ( fSplittingFactor  );
  fSplitAndKillOperation->SetApplyProbability( fApplyProbability );
  G4cout << GetName() << " : starting run with splitting factor = " << fSplittingFactor
         << ", and probability for applying the technique " << fApplyProbability
         << " . " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation*
GB03BOptrGeometryBasedBiasing::
ProposeNonPhysicsBiasingOperation( const G4Track*                   /* track */,
                                   const G4BiasingProcessInterface* /* callingProcess */ )
{
  // Here, we always return the split and kill biasing operation:
  return fSplitAndKillOperation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
