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
/// \file GB04BOptrBremSplitting.cc
/// \brief Implementation of the GB04BOptrBremSplitting class

#include "GB04BOptrBremSplitting.hh"
#include "GB04BOptnBremSplitting.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB04BOptrBremSplitting::GB04BOptrBremSplitting()
: G4VBiasingOperator("BremSplittingOperator"),
  fSplittingFactor(1),
  fBiasPrimaryOnly(true),
  fBiasOnlyOnce(true)
{
  fBremSplittingOperation = new GB04BOptnBremSplitting("BremSplittingOperation");
  
  // -- Define messengers:
  // -- Splitting factor:
  fSplittingFactorMessenger = 
    new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
  G4GenericMessenger::Command& splittingFactorCmd = 
    fSplittingFactorMessenger->DeclareProperty("setSplittingFactor", fSplittingFactor,
                                               "Define the brem. splitting factor." );
  splittingFactorCmd.SetStates(G4State_Idle);
  // -- Bias ony primary particle:
  fBiasPrimaryOnlyMessenger = 
    new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
  G4GenericMessenger::Command& biasPrimaryCmd = 
    fBiasPrimaryOnlyMessenger->DeclareProperty("biasPrimaryOnly", fBiasPrimaryOnly,
                      "Chose if brem. splitting applies to primary particles only." );
  biasPrimaryCmd.SetStates(G4State_Idle);
  // -- Bias ony primary particle:
  fBiasOnlyOnceMessenger = 
    new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
  G4GenericMessenger::Command& biasOnlyOnceCmd = 
    fBiasPrimaryOnlyMessenger->DeclareProperty("biasOnlyOnce", fBiasOnlyOnce,
                      "Chose if apply the brem. splitting only once for the track." );
  biasOnlyOnceCmd.SetStates(G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB04BOptrBremSplitting::StartRun()
{
  fBremSplittingOperation->SetSplittingFactor ( fSplittingFactor );
  G4cout << GetName() << " : starting run with brem. splitting factor = " 
         << fSplittingFactor;
  if ( fBiasPrimaryOnly ) G4cout << ", biasing only primaries ";
  else                    G4cout << ", biasing primary and secondary tracks ";
  if ( fBiasOnlyOnce )    G4cout << ", biasing only once per track ";
  else                    G4cout << ", biasing several times per track ";
  G4cout << " . " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB04BOptrBremSplitting::StartTracking( const G4Track* /* track */ )
{
  // -- reset the number of times the brem. splitting was applied:
  fNInteractions = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation*
GB04BOptrBremSplitting::
ProposeFinalStateBiasingOperation(const G4Track* track,
                                  const G4BiasingProcessInterface* /* callingProcess */)
{
  // -- Check if biasing of primary particle only is requested. If so, and
  // -- if particle is not a primary one, don't ask for biasing:
  if ( fBiasPrimaryOnly && ( track->GetParentID() !=0 ) ) return 0;
  // -- Check if brem. splitting should be applied only once to the track,
  // -- and if so, and if brem. splitting already occured, don't ask for biasing:
  if ( fBiasOnlyOnce    && ( fNInteractions > 0 )        ) return 0;
  
  // -- Count the number of times the brem. splitting is applied:
  fNInteractions++;
  // -- Return the brem. splitting operation:
  return fBremSplittingOperation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
