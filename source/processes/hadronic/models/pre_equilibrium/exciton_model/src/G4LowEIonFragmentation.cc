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
// ClassName:   G4LowEIonFragmentation
//
// Author:  H.P. Wellisch
//
// Modified:
// 02 Jun 2010 M. A. Cortes Giraldo fix: particlesFromTarget must be 
//                     accounted for as particles of initial compound nucleus
// 28 Oct 2010 V.Ivanchenko complete migration to integer Z and A; 
//                          use updated G4Fragment methods

#include <algorithm>

#include "G4LowEIonFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4Proton.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicsModelCatalog.hh"

G4LowEIonFragmentation::G4LowEIonFragmentation(G4ExcitationHandler * const value) 
  : G4HadronicInteraction("LowEIonPreco")
{
  theHandler = value;
  theModel = new G4PreCompoundModel(theHandler);
  proton = G4Proton::Proton();
  secID = G4PhysicsModelCatalog::GetModelID("model_" + GetModelName());
}

G4LowEIonFragmentation::~G4LowEIonFragmentation() 
{
  theResult.Clear();
}

G4HadFinalState * G4LowEIonFragmentation::
ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus)
{
  area = 0.0;
  // initialize the particle change
  theResult.Clear();
  theResult.SetStatusChange( stopAndKill );
  theResult.SetEnergyChange( 0.0 );

  // Get Target A, Z
  G4int aTargetA = theNucleus.GetA_asInt();
  G4int aTargetZ = theNucleus.GetZ_asInt();

  // Get Projectile A, Z
  G4int aProjectileA = thePrimary.GetDefinition()->GetBaryonNumber();
  G4int aProjectileZ = 
    G4lrint(thePrimary.GetDefinition()->GetPDGCharge()/eplus);

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

  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while(0==particlesFromProjectile)
  {
    do
    {
      x = 2*G4UniformRand() - 1;
      y = 2*G4UniformRand() - 1;
    }
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    while(x*x + y*y > 1);
    impactParameter = std::sqrt(x*x+y*y)*
      (targetOuterRadius+projectileOuterRadius);
    ++totalTries;
    area = pi*(targetOuterRadius+projectileOuterRadius)*
              (targetOuterRadius+projectileOuterRadius);
    G4double projectileHorizon = impactParameter-targetOuterRadius; 
    
    // Empirical boundary transparency.
    G4double empirical = G4UniformRand();
    if(projectileHorizon > empirical*projectileOuterRadius) { continue; }
    
    // Calculate the number of nucleons involved in collision
    // From projectile
    aPrim.StartLoop();

    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    while((pNucleon = aPrim.GetNextNucleon()))
    {
      if(pNucleon->GetPosition().y()>projectileHorizon)
      {
        // We have one
        ++particlesFromProjectile;
        if(pNucleon->GetParticleType() == proton) 
        {
          ++chargedFromProjectile;
        } 
      }
    }
  }
  ++hits;

  // From target:
  G4double targetHorizon = impactParameter-projectileOuterRadius;
  G4int chargedFromTarget = 0;
  G4int particlesFromTarget = 0;
  aTarg.StartLoop();  
  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while((pNucleon = aTarg.GetNextNucleon()))
  {
    if(pNucleon->GetPosition().y()>targetHorizon)
    {
      // We have one
      ++particlesFromTarget;
      if(pNucleon->GetParticleType() == proton) 
      {
        ++chargedFromTarget;
      }
    }
  }
  
  // Energy sharing between projectile and target. 
  // Note that this is a quite simplistic kinetically.
  G4ThreeVector momentum = thePrimary.Get4Momentum().vect();
  G4double w = (G4double)particlesFromProjectile/(G4double)aProjectileA;
  
  G4double projTotEnergy = thePrimary.GetTotalEnergy();  
  G4double targetMass = G4NucleiProperties::GetNuclearMass(aTargetA, aTargetZ);
  G4LorentzVector fragment4Momentum(momentum*w, projTotEnergy*w + targetMass);
 
  // take the nucleons and fill the Fragments
  G4Fragment anInitialState(aTargetA+particlesFromProjectile,
			    aTargetZ+chargedFromProjectile,
  			    fragment4Momentum);
  // M.A. Cortes fix
  anInitialState.SetNumberOfExcitedParticle(particlesFromProjectile
					    + particlesFromTarget,
					    chargedFromProjectile 
					    + chargedFromTarget);
  anInitialState.SetNumberOfHoles(particlesFromProjectile+particlesFromTarget,
				  chargedFromProjectile + chargedFromTarget);
  G4double time = thePrimary.GetGlobalTime();
  anInitialState.SetCreationTime(time);
  anInitialState.SetCreatorModelID(secID);

  // Fragment the Fragment using Pre-compound
  G4ReactionProductVector* thePreCompoundResult = 
    theModel->DeExcite(anInitialState);
  
  // De-excite the projectile using ExcitationHandler
  G4ReactionProductVector * theExcitationResult = nullptr;
  if(particlesFromProjectile < aProjectileA)
  {
    G4LorentzVector residual4Momentum(momentum*(1.0-w), projTotEnergy*(1.0-w));
 
    G4Fragment initialState2(aProjectileA-particlesFromProjectile,
			     aProjectileZ-chargedFromProjectile,
			     residual4Momentum );

    // half of particles are excited (?!)
    G4int pinit = (aProjectileA-particlesFromProjectile)/2;
    G4int cinit = (aProjectileZ-chargedFromProjectile)/2;

    initialState2.SetNumberOfExcitedParticle(pinit,cinit);
    initialState2.SetNumberOfHoles(pinit,cinit);
    initialState2.SetCreationTime(time);
    initialState2.SetCreatorModelID(secID);

    theExcitationResult = theHandler->BreakItUp(initialState2);
  }

  // Fill the particle change and clear intermediate vectors
  std::size_t nexc = (nullptr != theExcitationResult) ? 
    theExcitationResult->size() : 0;
  std::size_t npre = (nullptr != thePreCompoundResult) ?
    thePreCompoundResult->size() : 0;
  
  for(std::size_t k=0; k<nexc; ++k) {
    G4ReactionProduct* p = (*theExcitationResult)[k];
    G4HadSecondary secondary(new G4DynamicParticle(p->GetDefinition(), p->GetMomentum()));
    secondary.SetTime(p->GetTOF());
    secondary.SetCreatorModelID(secID);
    theResult.AddSecondary(secondary);
    delete p;
  }
  for(std::size_t k=0; k<npre; ++k) {
    G4ReactionProduct* p = (*thePreCompoundResult)[k];
    G4HadSecondary secondary(new G4DynamicParticle(p->GetDefinition(), p->GetMomentum()));
    secondary.SetTime(p->GetTOF());
    secondary.SetCreatorModelID(secID);
    theResult.AddSecondary(secondary);
    delete p;
  }
  
  delete thePreCompoundResult;
  delete theExcitationResult;

  // return the particle change
  return &theResult;
}
