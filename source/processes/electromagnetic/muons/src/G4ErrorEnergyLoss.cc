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

#include "G4ErrorEnergyLoss.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4ErrorPropagatorData.hh"

//-------------------------------------------------------------------
G4ErrorEnergyLoss::G4ErrorEnergyLoss(const G4String& processName, 
				     G4ProcessType type)
  : G4VContinuousProcess(processName, type)
{
  if (verboseLevel>2) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  theELossForExtrapolator = new G4EnergyLossForExtrapolator;
  theStepLimit = 1.*CLHEP::mm;
}

//-------------------------------------------------------------------
void G4ErrorEnergyLoss::InstantiateEforExtrapolator()
{}

//-------------------------------------------------------------------
G4ErrorEnergyLoss::~G4ErrorEnergyLoss() 
{
  delete theELossForExtrapolator;
}

//-------------------------------------------------------------------

G4bool G4ErrorEnergyLoss::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (aParticleType.GetPDGCharge() != 0);
}

//-------------------------------------------------------------------
G4VParticleChange*
G4ErrorEnergyLoss::AlongStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  G4ErrorPropagatorData* g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();

  G4double kinEnergyStart = aTrack.GetKineticEnergy();
  G4double step_length  = aStep.GetStepLength();

  const G4Material* aMaterial = aTrack.GetMaterial();
  const G4ParticleDefinition* aParticleDef = aTrack.GetDynamicParticle()->GetDefinition();
  G4double kinEnergyEnd = kinEnergyStart;

  // backward - energy increased
  if( g4edata->GetMode() == G4ErrorMode(G4ErrorMode_PropBackwards) ) {
    kinEnergyEnd = theELossForExtrapolator->EnergyBeforeStep( kinEnergyStart, 
							      step_length, 
							      aMaterial, 
							      aParticleDef );
    G4double kinEnergyHalfStep = (kinEnergyStart + kinEnergyEnd) * 0.5;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 ) 
    G4cout << " G4ErrorEnergyLoss FWD  end " << kinEnergyEnd 
	   << " halfstep " << kinEnergyHalfStep << G4endl;
#endif

    //--- rescale to energy lost at 1/2 step
    kinEnergyEnd = theELossForExtrapolator->EnergyBeforeStep( kinEnergyHalfStep, 
							      step_length, 
							      aMaterial, 
							      aParticleDef );
    kinEnergyEnd = kinEnergyStart - (kinEnergyHalfStep - kinEnergyEnd );

    // forward - energy decreased
  } else {
    kinEnergyEnd = theELossForExtrapolator->EnergyAfterStep( kinEnergyStart, 
							     step_length, 
							     aMaterial, 
							     aParticleDef );
    G4double kinEnergyHalfStep = (kinEnergyStart + kinEnergyEnd) * 0.5;
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 ) 
    G4cout << " G4ErrorEnergyLoss BCKD  end " << kinEnergyEnd 
	   << " halfstep " << kinEnergyHalfStep << G4endl;
#endif

    //--- rescale to energy lost at 1/2 step
    kinEnergyEnd = theELossForExtrapolator->EnergyAfterStep( kinEnergyHalfStep, 
							     step_length, 
							     aMaterial, 
							     aParticleDef );
    kinEnergyEnd = kinEnergyStart - (kinEnergyHalfStep - kinEnergyEnd );
  }

  G4double edepo = kinEnergyEnd - kinEnergyStart;

#ifdef G4VERBOSE
  if( G4ErrorPropagatorData::verbose() >= 2 ) 
    G4cout << "AlongStepDoIt Estart= " << kinEnergyStart << " Eend " << kinEnergyEnd 
	   << " Ediff " << kinEnergyStart-kinEnergyEnd << " step= " << step_length 
	   << " mate= " << aMaterial->GetName() 
	   << " particle= " << aParticleDef->GetParticleName() << G4endl;
#endif

  aParticleChange.ClearDebugFlag();
  aParticleChange.ProposeLocalEnergyDeposit( edepo );
  aParticleChange.SetNumberOfSecondaries(0);
 
  aParticleChange.ProposeEnergy( kinEnergyEnd );
  
  return &aParticleChange;
}


//-------------------------------------------------------------------
G4double G4ErrorEnergyLoss::GetContinuousStepLimit(const G4Track& aTrack,
				    G4double, G4double, G4double& )
{ 
  G4double ekin = aTrack.GetKineticEnergy();
  const G4Material* mat = aTrack.GetMaterial();
  const G4ParticleDefinition* part = 
    aTrack.GetDynamicParticle()->GetDefinition();
  G4double range = theELossForExtrapolator->ComputeRange(ekin, part, mat);
  G4double delta = std::max(range*theFractionLimit, theStepLimit);
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 ) {
    G4cout << " G4ErrorEnergyLoss: limiting Step " << delta 
	   << " energy(GeV) " << ekin / CLHEP::GeV 
	   << " for " << part->GetParticleName() << G4endl;
  }
#endif
  return delta;
}
