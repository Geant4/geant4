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
/// \file BiasingOperator.cc
/// \brief Implementation of the BiasingOperator class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BiasingOperator.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BiasingOperation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BiasingOperator::BiasingOperator() : G4VBiasingOperator( "BiasingOperator" ) {
  fBiasingOperation = new BiasingOperation( "BiasingOperation" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BiasingOperator::AddParticle( G4String particleName ) {
  const G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  if ( particle == 0 ) {
    G4ExceptionDescription ed;
    ed << "Particle `" << particleName << "' not found !" << G4endl;
    G4Exception( "BiasingOperator::AddParticle(...)", "BiasError", JustWarning, ed );
    return;
  }
  fParticlesToBias.push_back( particle );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation* BiasingOperator::
ProposeFinalStateBiasingOperation( const G4Track* , 
                                   const G4BiasingProcessInterface* callingProcess ) {
  // Apply the biasing operation only for inelastic processes of:
  // proton, neutron, pion+ and pion-
  if ( callingProcess  &&  callingProcess->GetWrappedProcess()  &&  
       ( callingProcess->GetWrappedProcess()->GetProcessName() == "protonInelastic"  ||
         callingProcess->GetWrappedProcess()->GetProcessName() == "neutronInelastic" || 
         callingProcess->GetWrappedProcess()->GetProcessName() == "pi+Inelastic"     || 
         callingProcess->GetWrappedProcess()->GetProcessName() == "pi-Inelastic" ) ) {
    return fBiasingOperation;
  } else {  
    return 0; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
