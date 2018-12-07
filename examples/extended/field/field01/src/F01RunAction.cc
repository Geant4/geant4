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
#include "F01RunAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"

#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

F01RunAction::F01RunAction() {
  theWarningEnergy   =   1.0 * CLHEP::kiloelectronvolt;  // Arbitrary 
  theImportantEnergy =  10.0 * CLHEP::kiloelectronvolt;  // Arbitrary 
  theNumberOfTrials  =   15;  // Arbitrary
  // Applications should determine these thresholds according to
  //  - physics requirements, and 
  //  - the computing cost of continuing integration for looping tracks
}

F01RunAction::~F01RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void F01RunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4cout << " Calling F01RunAction::ChangeLooperParameters() " << G4endl;
  ChangeLooperParameters( G4Electron::Definition() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void F01RunAction::
ChangeLooperParameters(const G4ParticleDefinition* particleDef )
{
  if( particleDef == nullptr )
     particleDef = G4Electron::Definition();
  auto transportPair= findTransportation( particleDef );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if( transport != nullptr )
  { 
     // Change the values of the looping particle parameters of Transportation 
     if( theWarningEnergy   >= 0.0 )
        transport->SetThresholdWarningEnergy(  theWarningEnergy ); 
     if( theImportantEnergy >= 0.0 )
        transport->SetThresholdImportantEnergy(  theImportantEnergy ); 
     if( theNumberOfTrials  > 0 )
        transport->SetThresholdTrials( theNumberOfTrials );
  }
  else if( coupledTransport != nullptr )
  { 
     // Change the values for Coupled Transport
     if( theWarningEnergy   >= 0.0 )
        coupledTransport->SetThresholdWarningEnergy(  theWarningEnergy ); 
     if( theImportantEnergy >= 0.0 )
        coupledTransport->SetThresholdImportantEnergy(  theImportantEnergy ); 
     if( theNumberOfTrials  > 0 )
        coupledTransport->SetThresholdTrials( theNumberOfTrials );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void F01RunAction::EndOfRunAction( const G4Run* ) {
  if( theVerboseLevel > 1 )
     G4cout << G4endl << G4endl
            << " ###########  Track Statistics for Transportation process(es) "
            << " ########### " << G4endl
            << " ############################################## "
            << " ####################### " << G4endl << G4endl;

  auto transportPair= findTransportation( G4Electron::Definition() );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;
  if( transport)             {        transport->PrintStatistics(G4cout); }
  else if( coupledTransport) { coupledTransport->PrintStatistics(G4cout); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::pair<G4Transportation*, G4CoupledTransportation*>
F01RunAction::findTransportation( const G4ParticleDefinition* particleDef,
                                  bool reportError )
{
  const auto *partPM=  particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport=
     dynamic_cast<G4CoupledTransportation*>(partTransport);

  if( reportError && !transport && !coupledTransport )
  {
     G4cerr << "Unable to find Transportation process for particle type "
            << particleDef->GetParticleName()
            << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
            << G4endl;
  }
  
  return
     std::make_pair( transport, coupledTransport );
         // <G4Transportation*, G4CoupledTransportation*>  
}
