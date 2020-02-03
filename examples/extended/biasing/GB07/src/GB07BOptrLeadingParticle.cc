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
/// \file GB07/src/GB07BOptrLeadingParticle.cc
/// \brief Implementation of the GB07BOptrLeadingParticle class
//
#include "GB07BOptrLeadingParticle.hh"
#include "G4BiasingProcessInterface.hh"

#include "G4BOptnLeadingParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PionZero.hh"
#include "G4ProcessManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07BOptrLeadingParticle::GB07BOptrLeadingParticle( G4String operatorName )
  : G4VBiasingOperator        ( operatorName ),
    fAnnihilation             ( nullptr ),
    fConversion               ( nullptr ),
    fDecay                    ( nullptr ),
    fTwoParticleProcess       ( nullptr )
{
  fLeadingParticleBiasingOperation =
    new G4BOptnLeadingParticle("LeadingParticleBiasingOperation");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07BOptrLeadingParticle::~GB07BOptrLeadingParticle()
{
  delete fLeadingParticleBiasingOperation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation* 
GB07BOptrLeadingParticle::
ProposeFinalStateBiasingOperation(const G4Track*, 
                                  const G4BiasingProcessInterface* callingProcess)
{
  // --   When the present method is called, we are at the process final state
  // -- generation level. The process is given by the callingProcess argument,
  // -- which, in our case,  wrappes a physics process, to control it.
  // --   To bias the final state generation, we return a biasing operation
  // -- which is fLeadingParticleBiasingOperation here. Before returning it, we
  // -- configure it depending on if the process is a two-particle final state
  // -- or if it is a many-particle final state process. For the two-particle
  // -- final state, one track is the leading and the other is alone in its category,
  // -- so always surviving by default. We play a Russian roulette on it to
  // -- trim also these two-particles final states.
  
  if ( callingProcess == fTwoParticleProcess )
    {
      // -- secondary particle accompagnying the leading one will be
      // -- killed with 2./3. probability (Russian roulette):
      fLeadingParticleBiasingOperation->SetFurtherKillingProbability( 2./3.);
    }
  else
    {
      // -- -1.0 means no effect : no further killing is applied to secondary
      // -- particles accompanying the leading one.
      fLeadingParticleBiasingOperation->SetFurtherKillingProbability( -1.0 );
    }
  
  return fLeadingParticleBiasingOperation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
GB07BOptrLeadingParticle::
StartRun()
{
  // -- collect the two-particle final state processes:
  fAnnihilation = nullptr;
  fConversion   = nullptr;
  fDecay        = nullptr;
  
  // ---- collect e+ annihilation process:
  auto positronProcesses = G4Positron::Definition()->GetProcessManager()->GetProcessList();
  for ( size_t i = 0; i < positronProcesses->size(); ++i )
    {
      if ( (*positronProcesses)[i]->GetProcessName() == "biasWrapper(annihil)")
        {
          fAnnihilation = (*positronProcesses)[i];
          break;
        }
    }
  
  // ---- collect gamma conversion process:
  auto gammaProcesses = G4Gamma::Definition()->GetProcessManager()->GetProcessList();
  for ( size_t i = 0; i < gammaProcesses->size(); ++i )
    {
      if ( (*gammaProcesses)[i]->GetProcessName() == "biasWrapper(conv)")
        {
          fConversion = (*gammaProcesses)[i];
          break;
        }
    }
  
  // ---- collect pi0 decay process:
  auto pi0Processes = G4PionZero::Definition()->GetProcessManager()->GetProcessList();
  for ( size_t i = 0; i < pi0Processes->size(); ++i )
    {
      if ( (*pi0Processes)[i]->GetProcessName() == "biasWrapper(Decay)")
        {
          fDecay = (*pi0Processes)[i];
          break;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
GB07BOptrLeadingParticle::
StartTracking( const G4Track* track )
{
  // -- remember what is the two-particle final state process -if any- for this starting
  // -- track:
  fTwoParticleProcess = nullptr;
  if ( track->GetDefinition() == G4Gamma   ::Definition() ) fTwoParticleProcess =   fConversion;
  if ( track->GetDefinition() == G4Positron::Definition() ) fTwoParticleProcess = fAnnihilation;
  if ( track->GetDefinition() == G4PionZero::Definition() ) fTwoParticleProcess =        fDecay;
}
