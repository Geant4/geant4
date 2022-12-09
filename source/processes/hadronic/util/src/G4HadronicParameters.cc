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
//---------------------------------------------------------------------------
//
// ClassName:      G4HadronicParameters
//
// Author:         2018 Alberto Ribon
//
// Description:    Singleton to keep global hadronic parameters.
//
// Modified:
//
//----------------------------------------------------------------------------

#include "G4HadronicParameters.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4HadronicParametersMessenger.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

G4HadronicParameters* G4HadronicParameters::sInstance = nullptr;

namespace
{
  G4Mutex paramMutex = G4MUTEX_INITIALIZER;
}

G4HadronicParameters* G4HadronicParameters::Instance() {
  if ( sInstance == nullptr ) {
    G4AutoLock l(&paramMutex);
    if ( sInstance == nullptr ) {
      static G4HadronicParameters theHadronicParametersObject;
      sInstance = &theHadronicParametersObject;
    }
    l.unlock();
  }
  return sInstance;
}


G4HadronicParameters::~G4HadronicParameters() {
  delete fMessenger;
}


G4HadronicParameters::G4HadronicParameters() {
  fMaxEnergy = 100.0*CLHEP::TeV;
  fMinEnergyTransitionFTF_Cascade = 3.0*CLHEP::GeV;
  fMaxEnergyTransitionFTF_Cascade = 6.0*CLHEP::GeV;
  fMinEnergyTransitionQGS_FTF = 12.0*CLHEP::GeV;
  fMaxEnergyTransitionQGS_FTF = 25.0*CLHEP::GeV;
  fEnergyThresholdForHeavyHadrons = 1.1*CLHEP::GeV;
  fMessenger = new G4HadronicParametersMessenger( this );
}


G4bool G4HadronicParameters::IsLocked() const {
  return ( ! G4Threading::IsMasterThread() ||
           G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit );
}


void G4HadronicParameters::SetMaxEnergy( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMaxEnergy = val;
  }
}


void G4HadronicParameters::SetMinEnergyTransitionFTF_Cascade( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMinEnergyTransitionFTF_Cascade = val;
  }
}


void G4HadronicParameters::SetMaxEnergyTransitionFTF_Cascade( const G4double val ) {
  if ( ! IsLocked()  &&  val > fMinEnergyTransitionFTF_Cascade ) { 
    fMaxEnergyTransitionFTF_Cascade = val;
  }
}


void G4HadronicParameters::SetMinEnergyTransitionQGS_FTF( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMinEnergyTransitionQGS_FTF = val;
  }
}

void G4HadronicParameters::SetMaxEnergyTransitionQGS_FTF( const G4double val ) {
  if ( ! IsLocked()  &&  val > fMinEnergyTransitionQGS_FTF ) { 
    fMaxEnergyTransitionQGS_FTF = val;
  }
}

void G4HadronicParameters::SetEnableBCParticles( G4bool val ) {
  if ( ! IsLocked() ) fEnableBC = val;
}


void G4HadronicParameters::SetEnableHyperNuclei( G4bool val ) {
  if ( ! IsLocked() ) fEnableHyperNuclei = val;
}


void G4HadronicParameters::SetVerboseLevel( const G4int val ) {
  if ( ! IsLocked()  &&  val >= 0 ) fVerboseLevel = val;
}


void G4HadronicParameters::SetEnergyThresholdForHeavyHadrons( G4double val ) {
  if ( ! IsLocked()  &&  val >= 0 && val < 5*CLHEP::GeV ) {
    fEnergyThresholdForHeavyHadrons = val;
  }
}


void G4HadronicParameters::SetXSFactorNucleonInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorNucleonInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorNucleonElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorNucleonElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorPionInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorPionInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorPionElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorPionElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorHadronInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorHadronInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorHadronElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorHadronElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorEM( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorEM = val;
  }
}


void G4HadronicParameters::SetApplyFactorXS( G4bool val ) {
  if ( ! IsLocked() ) fApplyFactorXS = val; 
}


void G4HadronicParameters::SetEnableCRCoalescence( G4bool val ) {
  if ( ! IsLocked() ) fEnableCRCoalescence = val;
}


void G4HadronicParameters::SetEnableIntegralInelasticXS( G4bool val ) {
  if ( ! IsLocked() ) fEnableIntegralInelasticXS = val;
}


void G4HadronicParameters::SetEnableIntegralElasticXS( G4bool val ) {
  if ( ! IsLocked() ) fEnableIntegralElasticXS = val;
}


void G4HadronicParameters::SetEnableDiffDissociationForBGreater10( G4bool val ) {
  if ( ! IsLocked() ) fEnableDiffDissociationForBGreater10 = val;
}


void G4HadronicParameters::SetEnableNeutronGeneralProcess( G4bool val ) {
  if ( ! IsLocked() ) fNeutronGeneral = val;
} 
