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
// $Id: G4RPGNeutronInelastic.cc 92494 2015-09-02 07:19:42Z gcosmo $
//
 
#include "G4RPGNeutronInelastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4HadFinalState* 
G4RPGNeutronInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                     G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile* originalIncident = &aTrack;

  // create the target particle
  G4DynamicParticle* originalTarget = targetNucleus.ReturnTargetParticle();
  
  G4ReactionProduct modifiedOriginal;
  modifiedOriginal = *originalIncident;
  G4ReactionProduct targetParticle;
  targetParticle = *originalTarget;
  if( originalIncident->GetKineticEnergy()/GeV < 0.01 + 2.*G4UniformRand()/9. )
  {
    SlowNeutron(originalIncident,modifiedOriginal,targetParticle,targetNucleus );
    delete originalTarget;
    return &theParticleChange;
  }

  // Fermi motion and evaporation
  // As of Geant3, the Fermi energy calculation had not been Done
  G4double ek = originalIncident->GetKineticEnergy()/MeV;
  G4double amas = originalIncident->GetDefinition()->GetPDGMass()/MeV;
    
  G4double tkin = targetNucleus.Cinema( ek );
  ek += tkin;
  modifiedOriginal.SetKineticEnergy( ek*MeV );
  G4double et = ek + amas;
  G4double p = std::sqrt( std::abs((et-amas)*(et+amas)) );
  G4double pp = modifiedOriginal.GetMomentum().mag()/MeV;
  if( pp > 0.0 )
  {
    G4ThreeVector momentum = modifiedOriginal.GetMomentum();
    modifiedOriginal.SetMomentum( momentum * (p/pp) );
  }
  //
  // calculate black track energies
  //
  tkin = targetNucleus.EvaporationEffects( ek );
  ek -= tkin;
  modifiedOriginal.SetKineticEnergy(ek);
  et = ek + amas;
  p = std::sqrt( std::abs((et-amas)*(et+amas)) );
  pp = modifiedOriginal.GetMomentum().mag();
  if( pp > 0.0 )
  {
    G4ThreeVector momentum = modifiedOriginal.GetMomentum();
    modifiedOriginal.SetMomentum( momentum * (p/pp) );
  }
  const G4double cutOff = 0.1;
  if( modifiedOriginal.GetKineticEnergy()/MeV <= cutOff )
  {
    SlowNeutron( originalIncident, modifiedOriginal, targetParticle, targetNucleus );
    delete originalTarget;
    return &theParticleChange;
  }

  G4ReactionProduct currentParticle = modifiedOriginal;
  currentParticle.SetSide( 1 ); // incident always goes in forward hemisphere
  targetParticle.SetSide( -1 );  // target always goes in backward hemisphere
  G4bool incidentHasChanged = false;
  G4bool targetHasChanged = false;
  G4bool quasiElastic = false;
  G4FastVector<G4ReactionProduct,256> vec;  // vec will contain sec. particles
  G4int vecLen = 0;
  vec.Initialize( 0 );
    
  InitialCollision(vec, vecLen, currentParticle, targetParticle,
                   incidentHasChanged, targetHasChanged);
    
  CalculateMomenta(vec, vecLen,
                   originalIncident, originalTarget, modifiedOriginal,
                   targetNucleus, currentParticle, targetParticle,
                   incidentHasChanged, targetHasChanged, quasiElastic);
    
  SetUpChange(vec, vecLen,
              currentParticle, targetParticle,
              incidentHasChanged);
    
  delete originalTarget;
  return &theParticleChange;
}
 
void
G4RPGNeutronInelastic::SlowNeutron(const G4HadProjectile* originalIncident,
                              G4ReactionProduct& modifiedOriginal,
                              G4ReactionProduct& targetParticle,
                              G4Nucleus& targetNucleus)
{        
  const G4double A = targetNucleus.GetA_asInt();    // atomic weight
  const G4double Z = targetNucleus.GetZ_asInt();    // atomic number
    
  G4double currentKinetic = modifiedOriginal.GetKineticEnergy()/MeV;
  G4double currentMass = modifiedOriginal.GetMass()/MeV;
  if( A < 1.5 )   // Hydrogen
  {
    //
    // very simple simulation of scattering angle and energy
    // nonrelativistic approximation with isotropic angular
    // distribution in the cms system 
    //
    G4double cost1, eka = 0.0;
    while (eka <= 0.0) {  /* Loop checking, 01.09.2015, D.Wright */
      cost1 = -1.0 + 2.0*G4UniformRand();
      eka = 1.0 + 2.0*cost1*A + A*A;
    }
    G4double cost = std::min( 1.0, std::max( -1.0, (A*cost1+1.0)/std::sqrt(eka) ) );
    eka /= (1.0+A)*(1.0+A);
    G4double ek = currentKinetic*MeV/GeV;
    G4double amas = currentMass*MeV/GeV;
    ek *= eka;
    G4double en = ek + amas;
    G4double p = std::sqrt(std::abs(en*en-amas*amas));
    G4double sint = std::sqrt(std::abs(1.0-cost*cost));
    G4double phi = G4UniformRand()*twopi;
    G4double px = sint*std::sin(phi);
    G4double py = sint*std::cos(phi);
    G4double pz = cost;
    targetParticle.SetMomentum( px*GeV, py*GeV, pz*GeV );
    G4double pxO = originalIncident->Get4Momentum().x()/GeV;
    G4double pyO = originalIncident->Get4Momentum().y()/GeV;
    G4double pzO = originalIncident->Get4Momentum().z()/GeV;
    G4double ptO = pxO*pxO + pyO+pyO;
    if( ptO > 0.0 )
    {
      G4double pO = std::sqrt(pxO*pxO+pyO*pyO+pzO*pzO);
      cost = pzO/pO;
      sint = 0.5*(std::sqrt(std::abs((1.0-cost)*(1.0+cost)))+std::sqrt(ptO)/pO);
      G4double ph = pi/2.0;
      if( pyO < 0.0 )ph = ph*1.5;
      if( std::abs(pxO) > 0.000001 )ph = std::atan2(pyO,pxO);
      G4double cosp = std::cos(ph);
      G4double sinp = std::sin(ph);
      px = cost*cosp*px - sinp*py+sint*cosp*pz;
      py = cost*sinp*px + cosp*py+sint*sinp*pz;
      pz = -sint*px     + cost*pz;
    }
    else
    {
      if( pz < 0.0 )pz *= -1.0;
    }
    G4double pu = std::sqrt(px*px+py*py+pz*pz);
    modifiedOriginal.SetMomentum( targetParticle.GetMomentum() * (p/pu) );
    modifiedOriginal.SetKineticEnergy( ek*GeV );
      
    targetParticle.SetMomentum(
     originalIncident->Get4Momentum().vect() - modifiedOriginal.GetMomentum() );
    G4double pp = targetParticle.GetMomentum().mag();
    G4double tarmas = targetParticle.GetMass();
    targetParticle.SetTotalEnergy( std::sqrt( pp*pp + tarmas*tarmas ) );
      
    theParticleChange.SetEnergyChange( modifiedOriginal.GetKineticEnergy() );
    G4DynamicParticle *pd = new G4DynamicParticle;
    pd->SetDefinition( targetParticle.GetDefinition() );
    pd->SetMomentum( targetParticle.GetMomentum() );
    theParticleChange.AddSecondary( pd );
    return;
  }

  G4FastVector<G4ReactionProduct,4> vec;  // vec will contain the secondary particles
  G4int vecLen = 0;
  vec.Initialize( 0 );
    
  G4double theAtomicMass = targetNucleus.AtomicMass( A, Z );
  G4double massVec[9];
  massVec[0] = targetNucleus.AtomicMass( A+1.0, Z     );
  massVec[1] = theAtomicMass;
  massVec[2] = 0.;
  if (Z > 1.0) massVec[2] = targetNucleus.AtomicMass(A, Z-1.0);
  massVec[3] = 0.;
  if (Z > 1.0 && A > 1.0) massVec[3] = targetNucleus.AtomicMass(A-1.0, Z-1.0 );
  massVec[4] = 0.;
  if (Z > 1.0 && A > 2.0 && A-2.0 > Z-1.0)
    massVec[4] = targetNucleus.AtomicMass( A-2.0, Z-1.0 );
  massVec[5] = 0.;
  if (Z > 2.0 && A > 3.0 && A-3.0 > Z-2.0)
    massVec[5] = targetNucleus.AtomicMass( A-3.0, Z-2.0 );
  massVec[6] = 0.;
  if (A > 1.0 && A-1.0 > Z) massVec[6] = targetNucleus.AtomicMass(A-1.0, Z);
  massVec[7] = massVec[3];
  massVec[8] = 0.;
  if (Z > 2.0 && A > 1.0) massVec[8] = targetNucleus.AtomicMass( A-1.0,Z-2.0 );
    
  twoBody.NuclearReaction(vec, vecLen, originalIncident,
                          targetNucleus, theAtomicMass, massVec );
    
  theParticleChange.SetStatusChange( stopAndKill );
  theParticleChange.SetEnergyChange( 0.0 );
    
  G4DynamicParticle* pd;
  for( G4int i=0; i<vecLen; ++i ) {
    pd = new G4DynamicParticle();
    pd->SetDefinition( vec[i]->GetDefinition() );
    pd->SetMomentum( vec[i]->GetMomentum() );
    theParticleChange.AddSecondary( pd );
    delete vec[i];
  }
}


// Initial Collision
//   selects the particle types arising from the initial collision of
//   the neutron and target nucleon.  Secondaries are assigned to
//   forward and backward reaction hemispheres, but final state energies
//   and momenta are not calculated here.

void
G4RPGNeutronInelastic::InitialCollision(G4FastVector<G4ReactionProduct,256>& vec,
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
  G4int part1;
  G4int part2;
 
  G4double testCharge;
  G4double testBaryon;
  G4double testStrange;
  
  // Get particle types according to incident and target types

  if (targetParticle.GetDefinition() == particleDef[neu]) {
    mult = GetMultiplicityT1(KE);
    fsTypes = GetFSPartTypesForNN(mult, KE);

    part1 = fsTypes[0];
    part2 = fsTypes[1];
    currentParticle.SetDefinition(particleDef[part1]);
    targetParticle.SetDefinition(particleDef[part2]);
    if (part1 == pro) {      
      if (part2 == neu) {
        if (G4UniformRand() > 0.5) {
          incidentHasChanged = true;
	} else {
          targetHasChanged = true;
          currentParticle.SetDefinition(particleDef[part2]);
          targetParticle.SetDefinition(particleDef[part1]);
	}
      } else {
        targetHasChanged = true;
        incidentHasChanged = true;
      }

    } else { // neutron
      if (part2 > neu && part2 < xi0) targetHasChanged = true;
    }

    testCharge = 0.0;
    testBaryon = 2.0;
    testStrange = 0.0;

  } else { // target was a proton
    mult = GetMultiplicityT0(KE);
    fsTypes = GetFSPartTypesForNP(mult, KE);
 
    part1 = fsTypes[0];
    part2 = fsTypes[1];
    currentParticle.SetDefinition(particleDef[part1]);
    targetParticle.SetDefinition(particleDef[part2]);
    if (part1 == pro) {
      if (part2 == pro) {
        incidentHasChanged = true;
      } else if (part2 == neu) {
        if (G4UniformRand() > 0.5) {
          incidentHasChanged = true;
          targetHasChanged = true;
	} else {
          currentParticle.SetDefinition(particleDef[part2]);
          targetParticle.SetDefinition(particleDef[part1]);
	}
        
      } else if (part2 > neu && part2 < xi0) {
        incidentHasChanged = true;
        targetHasChanged = true;
      }

    } else { // neutron
      targetHasChanged = true;
    }

    testCharge = 1.0;
    testBaryon = 2.0;
    testStrange = 0.0;
  }

  //  if (mult == 2 && !incidentHasChanged && !targetHasChanged)
  //                                              quasiElastic = true;
  
  // Remove incident and target from fsTypes
  
  fsTypes.erase(fsTypes.begin());
  fsTypes.erase(fsTypes.begin());

  // Remaining particles are secondaries.  Put them into vec.
  
  G4ReactionProduct* rp(0);
  for(G4int i=0; i < mult-2; ++i ) {
    partType = fsTypes[i];
    rp = new G4ReactionProduct();
    rp->SetDefinition(particleDef[partType]);
    (G4UniformRand() < 0.5) ? rp->SetSide(-1) : rp->SetSide(1);
    vec.SetElement(vecLen++, rp);
  }

  // Check conservation of charge, strangeness, baryon number
  
  CheckQnums(vec, vecLen, currentParticle, targetParticle,
             testCharge, testBaryon, testStrange);
  
  return;
}
