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
// $Id: G4RPGTwoBody.cc 94214 2015-11-09 08:18:05Z gcosmo $
//

#include <iostream>
#include <signal.h>

#include "G4RPGTwoBody.hh"
#include "G4Log.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4HadReentrentException.hh"

G4RPGTwoBody::G4RPGTwoBody()
  : G4RPGReaction() {}


G4bool G4RPGTwoBody::
ReactionStage(const G4HadProjectile* /*originalIncident*/,
              G4ReactionProduct& modifiedOriginal,
              G4bool& /*incidentHasChanged*/,
              const G4DynamicParticle* originalTarget,
              G4ReactionProduct& targetParticle,
              G4bool& /*targetHasChanged*/,
              const G4Nucleus& targetNucleus,
              G4ReactionProduct& currentParticle,
              G4FastVector<G4ReactionProduct,256>& vec,
              G4int& vecLen,
              G4bool /*leadFlag*/,
              G4ReactionProduct& /*leadingStrangeParticle*/)
{
  // 
  // Derived from H. Fesefeldt's original FORTRAN code TWOB
  //
  // Generation of momenta for elastic and quasi-elastic 2 body reactions
  //
  // The simple formula ds/d|t| = s0* std::exp(-b*|t|) is used.
  // The b values are parametrizations from experimental data.
  // Unavailable values are taken from those of similar reactions.
  //
    
  const G4double ekOriginal = modifiedOriginal.GetKineticEnergy()/GeV;
  const G4double etOriginal = modifiedOriginal.GetTotalEnergy()/GeV;
  const G4double mOriginal = modifiedOriginal.GetMass()/GeV;
  const G4double pOriginal = modifiedOriginal.GetMomentum().mag()/GeV;
  G4double currentMass = currentParticle.GetMass()/GeV;
  G4double targetMass = targetParticle.GetDefinition()->GetPDGMass()/GeV;

  targetMass = targetParticle.GetMass()/GeV;
  const G4double atomicWeight = targetNucleus.GetA_asInt();
    
  G4double etCurrent = currentParticle.GetTotalEnergy()/GeV;
  G4double pCurrent = currentParticle.GetTotalMomentum()/GeV;

  G4double cmEnergy = std::sqrt( currentMass*currentMass +
                            targetMass*targetMass +
                            2.0*targetMass*etCurrent );  // in GeV

  if (cmEnergy < 0.01) { // 2-body scattering not possible
    targetParticle.SetMass( 0.0 );  // flag that the target particle doesn't exist

  } else {
    // Projectile momentum in cm

    G4double pf = targetMass*pCurrent/cmEnergy;

    //
    // Set beam and target in centre of mass system
    //
    G4ReactionProduct pseudoParticle[3];
      
    if (targetParticle.GetDefinition()->GetParticleSubType() == "kaon" ||
        targetParticle.GetDefinition()->GetParticleSubType() == "pi") {

      //      G4double pf1 = pOriginal*mOriginal/std::sqrt(2.*mOriginal*(mOriginal+etOriginal));

      pseudoParticle[0].SetMass( targetMass*GeV );
      pseudoParticle[0].SetTotalEnergy( etOriginal*GeV );
      pseudoParticle[0].SetMomentum( 0.0, 0.0, pOriginal*GeV );
      
      pseudoParticle[1].SetMomentum( 0.0, 0.0, 0.0 );
      pseudoParticle[1].SetMass( mOriginal*GeV );
      pseudoParticle[1].SetKineticEnergy( 0.0 );

    } else {
      pseudoParticle[0].SetMass( currentMass*GeV );
      pseudoParticle[0].SetTotalEnergy( etCurrent*GeV );
      pseudoParticle[0].SetMomentum( 0.0, 0.0, pCurrent*GeV );
      
      pseudoParticle[1].SetMomentum( 0.0, 0.0, 0.0 );
      pseudoParticle[1].SetMass( targetMass*GeV );
      pseudoParticle[1].SetKineticEnergy( 0.0 );
    }
    //
    // Transform into center of mass system
    //
    pseudoParticle[2] = pseudoParticle[0] + pseudoParticle[1];
    pseudoParticle[0].Lorentz( pseudoParticle[0], pseudoParticle[2] );
    pseudoParticle[1].Lorentz( pseudoParticle[1], pseudoParticle[2] );
    //
    // Set final state masses and energies in centre of mass system
    //
    currentParticle.SetTotalEnergy( std::sqrt(pf*pf+currentMass*currentMass)*GeV );
    targetParticle.SetTotalEnergy( std::sqrt(pf*pf+targetMass*targetMass)*GeV );

    //
    // Calculate slope b for elastic scattering on proton/neutron
    //
    const G4double cb = 0.01;
    const G4double b1 = 4.225;
    const G4double b2 = 1.795;
    G4double b = std::max( cb, b1+b2*G4Log(pOriginal) );
     
    //
    // Get cm scattering angle by sampling t from tmin to tmax
    //
    G4double btrang = b * 4.0 * pf * pseudoParticle[0].GetMomentum().mag()/GeV;
    G4double exindt = G4Exp(-btrang) - 1.0;
    G4double costheta = 1.0 + 2*G4Log( 1.0+G4UniformRand()*exindt ) /btrang;
    costheta = std::max(-1., std::min(1., costheta) );
    G4double sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    G4double phi = twopi * G4UniformRand();

    //
    // Calculate final state momenta in centre of mass system
    //
    if (targetParticle.GetDefinition()->GetParticleSubType() == "kaon" ||
        targetParticle.GetDefinition()->GetParticleSubType() == "pi") {

      currentParticle.SetMomentum( -pf*sintheta*std::cos(phi)*GeV,
                                   -pf*sintheta*std::sin(phi)*GeV,
                                   -pf*costheta*GeV );
    } else {

      currentParticle.SetMomentum( pf*sintheta*std::cos(phi)*GeV,
                                   pf*sintheta*std::sin(phi)*GeV,
                                   pf*costheta*GeV );
    }
    targetParticle.SetMomentum( -currentParticle.GetMomentum() );

    //
    // Transform into lab system
    //
    currentParticle.Lorentz( currentParticle, pseudoParticle[1] );
    targetParticle.Lorentz( targetParticle, pseudoParticle[1] );

    // Rotate final state particle vectors wrt incident momentum
      
    Defs1( modifiedOriginal, currentParticle, targetParticle, vec, vecLen );

    // Subtract binding energy

    G4double pp, pp1, ekin;
    if( atomicWeight >= 1.5 )
    {
      const G4double cfa = 0.025*((atomicWeight-1.)/120.)*G4Exp(-(atomicWeight-1.)/120.);
      pp1 = currentParticle.GetMomentum().mag()/MeV;
      if( pp1 >= 1.0 )
      {
        ekin = currentParticle.GetKineticEnergy()/MeV - cfa*(1.0+0.5*normal())*GeV;
        ekin = std::max( 0.0001*GeV, ekin );
        currentParticle.SetKineticEnergy( ekin*MeV );
        pp = currentParticle.GetTotalMomentum()/MeV;
        currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
      }
      pp1 = targetParticle.GetMomentum().mag()/MeV;
      if( pp1 >= 1.0 )
      {
        ekin = targetParticle.GetKineticEnergy()/MeV - cfa*(1.0+normal()/2.)*GeV;
        ekin = std::max( 0.0001*GeV, ekin );
        targetParticle.SetKineticEnergy( ekin*MeV );
        pp = targetParticle.GetTotalMomentum()/MeV;
        targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
      }
    }
  }

  // Get number of final state nucleons and nucleons remaining in
  // target nucleus

  std::pair<G4int, G4int> finalStateNucleons = 
    GetFinalStateNucleons(originalTarget, vec, vecLen);
  G4int protonsInFinalState = finalStateNucleons.first;
  G4int neutronsInFinalState = finalStateNucleons.second;

  G4int PinNucleus = std::max(0, 
    G4int(targetNucleus.GetZ_asInt()) - protonsInFinalState);
  G4int NinNucleus = std::max(0,
    G4int(targetNucleus.GetA_asInt()-targetNucleus.GetZ_asInt()) - neutronsInFinalState);

  if( atomicWeight >= 1.5 )
  {
    // Add black track particles
    // npnb: number of proton/neutron black track particles
    // ndta: number of deuterons, tritons, and alphas produced
    // epnb: kinetic energy available for proton/neutron black track 
    //   particles
    // edta: kinetic energy available for deuteron/triton/alpha particles

    G4double epnb, edta;
    G4int npnb=0, ndta=0;
      
    epnb = targetNucleus.GetPNBlackTrackEnergy();   // was enp1 in fortran code
    edta = targetNucleus.GetDTABlackTrackEnergy();  // was enp3 in fortran code
    const G4double pnCutOff = 0.0001;       // GeV
    const G4double dtaCutOff = 0.0001;      // GeV
    //    const G4double kineticMinimum = 0.0001;
    //    const G4double kineticFactor = -0.010;
    //    G4double sprob = 0.0; // sprob = probability of self-absorption in heavy molecules
    if( epnb >= pnCutOff )
    {
      npnb = G4Poisson( epnb/0.02 );
      if( npnb > atomicWeight )npnb = G4int(atomicWeight);
      if( (epnb > pnCutOff) && (npnb <= 0) )npnb = 1;
      npnb = std::min( npnb, 127-vecLen );
    }
    if( edta >= dtaCutOff )
    {
      ndta = G4int(2.0 * G4Log(atomicWeight));
      ndta = std::min( ndta, 127-vecLen );
    }

    if (npnb == 0 && ndta == 0) npnb = 1;

    AddBlackTrackParticles(epnb, npnb, edta, ndta, modifiedOriginal, 
                           PinNucleus, NinNucleus, targetNucleus,
                           vec, vecLen);
  }

  //
  //  calculate time delay for nuclear reactions
  //
  if( (atomicWeight >= 1.5) && (atomicWeight <= 230.0) && (ekOriginal <= 0.2) )
    currentParticle.SetTOF( 1.0-500.0*G4Exp(-ekOriginal/0.04)*G4Log(G4UniformRand()) );
  else
    currentParticle.SetTOF( 1.0 );

  return true;
}
 
/* end of file */
