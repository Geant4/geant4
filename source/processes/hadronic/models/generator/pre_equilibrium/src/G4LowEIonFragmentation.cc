//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4LowEIonFragmentation.hh"
#include "g4std/algorithm"

G4int G4LowEIonFragmentation::hits = 0;
G4int G4LowEIonFragmentation::totalTries = 0;
G4double G4LowEIonFragmentation::area = 0;

G4VParticleChange * G4LowEIonFragmentation::
ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus)
{
  area = 0;
  // initialize the particle change
//  theResult.Initialize(thePrimary);
  theResult.SetStatusChange( fStopAndKill );
  theResult.SetEnergyChange( 0.0 );

  // Get Target A, Z
  G4double aTargetA = theNucleus.GetN();
  G4double aTargetZ = theNucleus.GetZ();

  // Get Projectile A, Z
  G4double aProjectileA = thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber();
  G4double aProjectileZ = thePrimary.GetDynamicParticle()->GetDefinition()->GetPDGCharge();

  // Get Maximum radius of both
  
  G4Fancy3DNucleus aPrim;
  aPrim.Init(aProjectileA, aProjectileZ);
  G4double projectileOuterRadius = aPrim.GetOuterRadius();
  
  G4Fancy3DNucleus aTarg;
  aTarg.Init(aTargetA, aTargetZ);
  G4double targetOuterRadius = aTarg.GetOuterRadius();

  // Get the Impact parameter
  G4int particlesFromProjectile = 0;
  G4int chargedFromProjectile = 0;
  G4double impactParameter = 0;
  G4double x,y;
  G4Nucleon * pNucleon;
  // need at lease one particle from the projectile model beyond the 
  // projectileHorizon.
  while(0==particlesFromProjectile)
  {
    do
    {
      x = 2*G4UniformRand() - 1;
      y = 2*G4UniformRand() - 1;
    }
    while(x*x + y*y > 1);
    impactParameter = sqrt(x*x+y*y)*(targetOuterRadius+projectileOuterRadius);
    totalTries++;
    area = pi*(targetOuterRadius+projectileOuterRadius)*
              (targetOuterRadius+projectileOuterRadius);
    G4double projectileHorizon = impactParameter-targetOuterRadius; 
    
    // Empirical boundary transparency.
    G4double empirical = G4UniformRand();
    if(projectileHorizon/projectileOuterRadius>empirical) continue;
    
    // Calculate the number of nucleons involved in collision
    // From projectile
    aPrim.StartLoop();
    while((pNucleon = aPrim.GetNextNucleon()))
    {
      if(pNucleon->GetPosition().y()>projectileHorizon)
      {
        // We have one
        particlesFromProjectile++;
        if(pNucleon->GetParticleType()==G4Proton::ProtonDefinition()) 
        {
          chargedFromProjectile++;
        } 
      }
    }
  }
  hits ++;

  // From target:
  G4double targetHorizon = impactParameter-projectileOuterRadius;
  G4int chargedFromTarget = 0;
  G4int particlesFromTarget = 0;
  aTarg.StartLoop();  
  while((pNucleon = aTarg.GetNextNucleon()))
  {
    if(pNucleon->GetPosition().y()>targetHorizon)
    {
      // We have one
      particlesFromTarget++;
      if(pNucleon->GetParticleType()==G4Proton::ProtonDefinition()) 
      {
        chargedFromTarget++;
      }
    }
  }
  
  // Energy sharing between projectile and target. Note that this is a quite simplistic kinetically.
  G4ThreeVector exciton3Momentum = thePrimary.GetMomentum();
  exciton3Momentum *= particlesFromProjectile/aProjectileA;
  
  G4double compoundEnergy = thePrimary.GetTotalEnergy()*particlesFromProjectile/aProjectileA;  
  G4double targetMass = G4ParticleTable::GetParticleTable()
                        ->GetIonTable()->GetIonMass(static_cast<G4int>(aTargetZ) ,static_cast<G4int>(aTargetA));
  compoundEnergy += targetMass;
  G4LorentzVector fragment4Momentum(exciton3Momentum, compoundEnergy);
 
  // take the nucleons and fill the Fragments
  G4Fragment anInitialState;
  anInitialState.SetA(aTargetA+particlesFromProjectile);
  anInitialState.SetZ(aTargetZ+chargedFromProjectile);
  anInitialState.SetNumberOfParticles(particlesFromProjectile);
  anInitialState.SetNumberOfHoles(particlesFromTarget);
  anInitialState.SetNumberOfCharged(chargedFromProjectile + chargedFromTarget);
  anInitialState.SetMomentum(fragment4Momentum);

  // Fragment the Fragment using Pre-compound
  G4ReactionProductVector* thePreCompoundResult;
  thePreCompoundResult = theModel->DeExcite(anInitialState);
  
  // De-excite the projectile using ExcitationHandler
  
  G4ReactionProductVector * theExcitationResult = 0; 
  if(particlesFromProjectile != aProjectileA)
  {
    G4ThreeVector residual3Momentum = thePrimary.GetMomentum();
    residual3Momentum -= exciton3Momentum;
    G4double residualEnergy = thePrimary.GetTotalEnergy()*(1.-particlesFromProjectile/aProjectileA);
    G4LorentzVector residual4Momentum(residual3Momentum, residualEnergy);  
 
    G4Fragment initialState2;
    initialState2.SetA(aProjectileA-particlesFromProjectile);
    initialState2.SetZ(aProjectileZ-chargedFromProjectile);
    initialState2.SetNumberOfHoles((aProjectileA-particlesFromProjectile)/2);
    initialState2.SetNumberOfParticles((aProjectileZ-chargedFromProjectile)/2);
    initialState2.SetNumberOfCharged((aProjectileZ-chargedFromProjectile)/2);


    initialState2.SetMomentum(residual4Momentum);
    theExcitationResult = theHandler->BreakItUp(initialState2);
  }

  // Fill the particle change
  G4int nSecondaries = 0;
  if(theExcitationResult) nSecondaries+=theExcitationResult->size();
  if(thePreCompoundResult) nSecondaries+=thePreCompoundResult->size();
  theResult.SetNumberOfSecondaries(nSecondaries);
  
  unsigned int k;
  if(theExcitationResult!=0)
  {
    for(k=0; k<theExcitationResult->size(); k++)
    {
      G4DynamicParticle* p0 = new G4DynamicParticle;
      p0->SetDefinition( theExcitationResult->operator[](k)->GetDefinition() );
      p0->SetMomentum( theExcitationResult->operator[](k)->GetMomentum() );
      theResult.AddSecondary(p0);
    }
  }
  
  for(k=0; k<thePreCompoundResult->size(); k++)
  {
    G4DynamicParticle* p0 = new G4DynamicParticle;
    p0->SetDefinition(thePreCompoundResult->operator[](k)->GetDefinition());
    p0->SetMomentum(thePreCompoundResult->operator[](k)->GetMomentum());
    theResult.AddSecondary(p0);
  }
  
  // clean up
  G4std::for_each(thePreCompoundResult->begin(), thePreCompoundResult->end(), DeleteReactionProduct());
  if(theExcitationResult) 
  {
    G4std::for_each(theExcitationResult->begin(), theExcitationResult->end(), DeleteReactionProduct());
  }
  delete thePreCompoundResult;
  if(theExcitationResult) delete theExcitationResult;

  // return the particle change
  return &theResult;
  
}
