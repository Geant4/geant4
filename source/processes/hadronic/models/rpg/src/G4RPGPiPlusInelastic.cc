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
// $Id: G4RPGPiPlusInelastic.cc 79697 2014-03-12 13:10:09Z gcosmo $
//
 
#include "G4RPGPiPlusInelastic.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4HadFinalState*
G4RPGPiPlusInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                     G4Nucleus& targetNucleus)
{
  const G4HadProjectile *originalIncident = &aTrack;
  if (originalIncident->GetKineticEnergy()<= 0.1) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit()); 
    return &theParticleChange;      
  }

    // create the target particle
    
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    G4ReactionProduct targetParticle( originalTarget->GetDefinition() );
    
    G4ReactionProduct currentParticle(originalIncident->GetDefinition() );
    currentParticle.SetMomentum( originalIncident->Get4Momentum().vect() );
    currentParticle.SetKineticEnergy( originalIncident->GetKineticEnergy() );
    
    // Fermi motion and evaporation
    // As of Geant3, the Fermi energy calculation had not been Done
    
    G4double ek = originalIncident->GetKineticEnergy();
    G4double amas = originalIncident->GetDefinition()->GetPDGMass();
    
    G4double tkin = targetNucleus.Cinema( ek );
    ek += tkin;
    currentParticle.SetKineticEnergy( ek );
    G4double et = ek + amas;
    G4double p = std::sqrt( std::abs((et-amas)*(et+amas)) );
    G4double pp = currentParticle.GetMomentum().mag();
    if( pp > 0.0 )
    {
      G4ThreeVector momentum = currentParticle.GetMomentum();
      currentParticle.SetMomentum( momentum * (p/pp) );
    }
    
    // calculate black track energies
    
    tkin = targetNucleus.EvaporationEffects( ek );
    ek -= tkin;
    currentParticle.SetKineticEnergy( ek );
    et = ek + amas;
    p = std::sqrt( std::abs((et-amas)*(et+amas)) );
    pp = currentParticle.GetMomentum().mag();
    if( pp > 0.0 )
    {
      G4ThreeVector momentum = currentParticle.GetMomentum();
      currentParticle.SetMomentum( momentum * (p/pp) );
    }

    G4ReactionProduct modifiedOriginal = currentParticle;

    currentParticle.SetSide( 1 ); // incident always goes in forward hemisphere
    targetParticle.SetSide( -1 );  // target always goes in backward hemisphere
    G4bool incidentHasChanged = false;
    G4bool targetHasChanged = false;
    G4bool quasiElastic = false;
    G4FastVector<G4ReactionProduct,256> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    const G4double cutOff = 0.1;
    if( currentParticle.GetKineticEnergy() > cutOff )
      InitialCollision(vec, vecLen, currentParticle, targetParticle,
                       incidentHasChanged, targetHasChanged);
    
    CalculateMomenta( vec, vecLen,
                      originalIncident, originalTarget, modifiedOriginal,
                      targetNucleus, currentParticle, targetParticle,
                      incidentHasChanged, targetHasChanged, quasiElastic );
    
    SetUpChange( vec, vecLen,
                 currentParticle, targetParticle,
                 incidentHasChanged );
    
    delete originalTarget;
    return &theParticleChange;
}


// Initial Collision
//   selects the particle types arising from the initial collision of
//   the projectile and target nucleon.  Secondaries are assigned to 
//   forward and backward reaction hemispheres, but final state energies
//   and momenta are not calculated here.
 
void 
G4RPGPiPlusInelastic::InitialCollision(G4FastVector<G4ReactionProduct,256>& vec,
                                  G4int& vecLen,
                                  G4ReactionProduct& currentParticle,
                                  G4ReactionProduct& targetParticle,
                                  G4bool& incidentHasChanged,
                                  G4bool& targetHasChanged)
{
  G4double KE = currentParticle.GetKineticEnergy()/GeV;

  G4int mult;
  G4int partType;
  std::vector<G4int> fsTypes;

  G4double testCharge;
  G4double testBaryon;
  G4double testStrange;

  // Get particle types according to incident and target types

  if (targetParticle.GetDefinition() == particleDef[pro]) {
    mult = GetMultiplicityT32(KE);
    fsTypes = GetFSPartTypesForPipP(mult, KE);
    partType = fsTypes[0];
    if (partType != pro) {
      targetHasChanged = true;
      targetParticle.SetDefinition(particleDef[partType]);
    }

    testCharge = 2.0;
    testBaryon = 1.0;
    testStrange = 0.0;

  } else {   // target was a neutron
    mult = GetMultiplicityT12(KE);
    fsTypes = GetFSPartTypesForPipN(mult, KE);
    partType = fsTypes[0];
    if (partType != neu) {
      targetHasChanged = true;
      targetParticle.SetDefinition(particleDef[partType]);
    }

    testCharge = 1.0;
    testBaryon = 1.0;
    testStrange = 0.0;
  }

  // Remove target particle from list

  fsTypes.erase(fsTypes.begin());

  // See if the incident particle changed type 

  G4int choose = -1;
  for(G4int i=0; i < mult-1; ++i ) {
    partType = fsTypes[i];
    if (partType == pip) {
      choose = i;
      break;
    }
  }
  if (choose == -1) {
    incidentHasChanged = true;
    choose = G4int(G4UniformRand()*(mult-1) );
    partType = fsTypes[choose];
    currentParticle.SetDefinition(particleDef[partType]);
  }
  fsTypes.erase(fsTypes.begin()+choose);

  // Remaining particles are secondaries.  Put them into vec.
  //   Improve this by randomizing secondary order, then alternate
  //   which secondary is put into forward or backward hemisphere

  G4ReactionProduct* rp(0);
  for(G4int i=0; i < mult-2; ++i ) {
    partType = fsTypes[i];
    rp = new G4ReactionProduct();
    rp->SetDefinition(particleDef[partType]);
    (G4UniformRand() < 0.5) ? rp->SetSide(-1) : rp->SetSide(1);
    if (partType > pim && partType < pro) rp->SetMayBeKilled(false);  // kaons
    vec.SetElement(vecLen++, rp);
  }
 
  //  if (mult == 2 && !incidentHasChanged && !targetHasChanged) 
  //                                              quasiElastic = true;

  // Check conservation of charge, strangeness, baryon number

  CheckQnums(vec, vecLen, currentParticle, targetParticle,
             testCharge, testBaryon, testStrange);

  return;
}
