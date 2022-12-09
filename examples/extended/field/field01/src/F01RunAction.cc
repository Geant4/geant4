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
  fWarningEnergy   =   1.0 * CLHEP::kiloelectronvolt;  // Arbitrary
  fImportantEnergy =  10.0 * CLHEP::kiloelectronvolt;  // Arbitrary
  fNumberOfTrials  =   15;  // Arbitrary
  // Applications should determine these thresholds according to
  //  - physics requirements, and 
  //  - the computing cost of continuing integration for looping tracks
}

F01RunAction::~F01RunAction() {}

//------------------------------------------------------------------------------


void F01RunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4cout << " Calling F01RunAction::ChangeLooperParameters() " << G4endl;
  ChangeLooperParameters( G4Electron::Definition() );
}

//------------------------------------------------------------------------------

void F01RunAction::
ChangeLooperParameters(const G4ParticleDefinition* particleDef )
{
  if( particleDef == nullptr )
     particleDef = G4Electron::Definition();
  auto transport= FindTransportation( particleDef );

  // Note that all particles that share the same Transportation
  //  (or Coupled Transportation) process (with this particle)
  //  will be affected by this change
  //   -- to change these parameters by particle type, the values
  //      must instead be changed/reset in a (Start) Tracking Action.

  G4cout << " ChangeLooperParameters called with particle type "
         << particleDef->GetParticleName()
         << " transport process: ";
  if( transport != nullptr ) {
    G4cout << transport->GetProcessName();
  } else {
    G4cout << " UNKNOWN -- it is neither G4Transportation nor G4CoupledTransportation";
  }
  G4cout << G4endl;
  
  if( transport != nullptr ) {
    if( fWarningEnergy   >= 0.0 ){
      transport->SetThresholdWarningEnergy(  fWarningEnergy );
      G4cout << "-- Changed Threshold Warning Energy (for loopers) = "
      << fWarningEnergy / CLHEP::MeV << " MeV "  << G4endl;
    }
    if( fImportantEnergy >= 0.0 ) {
      transport->SetThresholdImportantEnergy(  fImportantEnergy );
      
      G4cout << "-- Changed Threshold Important Energy (for loopers) = "
      << fImportantEnergy / CLHEP::MeV << " MeV " << G4endl;
    }
    
    if( fNumberOfTrials  > 0 ) {
      transport->SetThresholdTrials( fNumberOfTrials );
      
      G4cout << "-- Changed number of Trials (for loopers) = " << fNumberOfTrials << G4endl;
    }
  }
  
  if( transport == nullptr ) {
     if( fWarningEnergy   >= 0.0 )
        G4cerr << " Unknown transport process> Cannot change Warning Energy. " << G4endl;
     if( fImportantEnergy >= 0.0 )
        G4cerr << " Unknown transport process> Cannot change 'Important' Energy. " << G4endl;
     if( fNumberOfTrials  > 0 )
        G4cerr << " Unknown transport process> Cannot change number of trials. " << G4endl;
  }
}

//------------------------------------------------------------------------------

void F01RunAction::EndOfRunAction( const G4Run* ) {
  if( fVerboseLevel > 1 )
     G4cout << G4endl << G4endl
            << " ###########  Track Statistics for Transportation process(es) "
            << " ########### " << G4endl
            << " ############################################## "
            << " ####################### " << G4endl << G4endl;

  auto transport= FindTransportation( G4Electron::Definition() );
  if( transport) {
     transport->PrintStatistics(G4cout);
  }
}

//------------------------------------------------------------------------------

G4Transportation*
F01RunAction::FindTransportation( const G4ParticleDefinition* particleDef,
                                  bool reportError )
{
  const auto *partPM=  particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  if( !transport ){
     partTransport = partPM->GetProcess("CoupledTransportation");
     auto coupledTransport=
        dynamic_cast<G4CoupledTransportation*>(partTransport);
     if( coupledTransport ) {
        transport= coupledTransport;
     } else {
        partTransport = partPM->GetProcess("TransportationWithMsc");
        auto transportWithMsc=
           dynamic_cast<G4CoupledTransportation*>(partTransport);
        if( transportWithMsc ) {
           transport= transportWithMsc;
        }
     }
  }
  
  if( reportError && !transport )
  {
     G4cerr << "Unable to find Transportation process for particle type "
            << particleDef->GetParticleName()
            << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
            << G4endl;
  }
  
  return transport;
}
