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
// $Id: G4RPGFragmentation.cc 94214 2015-11-09 08:18:05Z gcosmo $
//
 
#include <iostream>
#include <signal.h>

#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4RPGFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4HadReentrentException.hh"

G4RPGFragmentation::G4RPGFragmentation()
 : G4RPGReaction()
{
  for (G4int i = 0; i < 20; i++) dndl[i] = 0.0;
}


void G4RPGFragmentation::
FragmentationIntegral(G4double pt, G4double et, G4double parMass, G4double secMass)
{
  pt = std::max( 0.001, pt );
  G4double dx = 1./(19.*pt);
  G4double x;
  G4double term1;
  G4double term2;

  for (G4int i = 1; i < 20; i++) {
    x = (G4double(i) - 0.5)*dx;
    term1 = 1. + parMass*parMass*x*x;
    term2 = pt*x*et*pt*x*et + pt*pt + secMass*secMass;
    dndl[i] = dx / std::sqrt( term1*term1*term1*term2 )
            + dndl[i-1];
  }
}


G4bool G4RPGFragmentation::
ReactionStage(const G4HadProjectile* originalIncident,
              G4ReactionProduct& modifiedOriginal,
              G4bool& incidentHasChanged,
              const G4DynamicParticle* originalTarget,
              G4ReactionProduct& targetParticle,
              G4bool& targetHasChanged,
              const G4Nucleus& targetNucleus,
              G4ReactionProduct& currentParticle,
              G4FastVector<G4ReactionProduct,256>& vec,
              G4int& vecLen,
              G4bool leadFlag,
              G4ReactionProduct& leadingStrangeParticle) 
{
  // 
  // Based on H. Fesefeldt's original FORTRAN code GENXPT
  //
  // Generation of x- and pT- values for incident, target, and all secondary 
  // particles using a simple single variable description E D3S/DP3= F(Q) 
  // with Q^2 = (M*X)^2 + PT^2.  Final state kinematics are produced by an 
  // FF-type iterative cascade method
  //
  // Internal units are GeV
  //
    
  // Protection in case no secondary has been created. In that case use  
  // two-body scattering
  //
  if (vecLen == 0) return false;

  G4ParticleDefinition* aPiMinus = G4PionMinus::PionMinus();
  G4ParticleDefinition* aProton = G4Proton::Proton();
  G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
  G4ParticleDefinition* aPiPlus = G4PionPlus::PionPlus();
  G4ParticleDefinition* aPiZero = G4PionZero::PionZero();

  G4int i, l;
  G4bool veryForward = false;
    
  const G4double ekOriginal = modifiedOriginal.GetKineticEnergy()/GeV;
  const G4double etOriginal = modifiedOriginal.GetTotalEnergy()/GeV;
  const G4double mOriginal = modifiedOriginal.GetMass()/GeV;
  const G4double pOriginal = modifiedOriginal.GetMomentum().mag()/GeV;
  G4double targetMass = targetParticle.GetDefinition()->GetPDGMass()/GeV;
  G4double centerofmassEnergy = std::sqrt( mOriginal*mOriginal +
                                      targetMass*targetMass +
                                      2.0*targetMass*etOriginal );  // GeV
  G4double currentMass = currentParticle.GetMass()/GeV;
  targetMass = targetParticle.GetMass()/GeV;

  // Randomize the order of the secondary particles.
  // Note that the current and target particles are not affected.

  for (i=0; i<vecLen; ++i) {
    G4int itemp = G4int( G4UniformRand()*vecLen );
    G4ReactionProduct pTemp = *vec[itemp];
    *vec[itemp] = *vec[i];
    *vec[i] = pTemp;
  }

  if (currentMass == 0.0 && targetMass == 0.0) {
    // Target and projectile have annihilated.  Replace them with the first 
    // two secondaries in the list.  Current particle KE is maintained.
 
    G4double ek = currentParticle.GetKineticEnergy();
    G4ThreeVector mom = currentParticle.GetMomentum();
    currentParticle = *vec[0];
    currentParticle.SetSide(1);
    targetParticle = *vec[1];
    targetParticle.SetSide(-1);
    for( i=0; i<(vecLen-2); ++i )*vec[i] = *vec[i+2];
    G4ReactionProduct *temp = vec[vecLen-1];
    delete temp;
    temp = vec[vecLen-2];
    delete temp;
    vecLen -= 2;
    currentMass = currentParticle.GetMass()/GeV;
    targetMass = targetParticle.GetMass()/GeV;
    incidentHasChanged = true;
    targetHasChanged = true;
    currentParticle.SetKineticEnergy( ek );
    currentParticle.SetMomentum(mom);
    veryForward = true;
  }
  const G4double atomicWeight = targetNucleus.GetA_asInt();
  const G4double atomicNumber = targetNucleus.GetZ_asInt();
  const G4double protonMass = aProton->GetPDGMass();

  if (originalIncident->GetDefinition()->GetParticleSubType() == "kaon"
      && G4UniformRand() >= 0.7) {
    G4ReactionProduct temp = currentParticle;
    currentParticle = targetParticle;
    targetParticle = temp;
    incidentHasChanged = true;
    targetHasChanged = true;
    currentMass = currentParticle.GetMass()/GeV;
    targetMass = targetParticle.GetMass()/GeV;
  }
  const G4double afc = std::min( 0.75,
        0.312+0.200*G4Log(G4Log(centerofmassEnergy*centerofmassEnergy))+
        std::pow(centerofmassEnergy*centerofmassEnergy,1.5)/6000.0 );
    
  G4double freeEnergy = centerofmassEnergy-currentMass-targetMass;
  G4double forwardEnergy = freeEnergy/2.;
  G4int forwardCount = 1;         // number of particles in forward hemisphere
    
  G4double backwardEnergy = freeEnergy/2.;
  G4int backwardCount = 1;        // number of particles in backward hemisphere

  if(veryForward) {
    if(currentParticle.GetSide()==-1)
    {
      forwardEnergy += currentMass;
      forwardCount --;
      backwardEnergy -= currentMass;
      backwardCount ++;
    }
    if(targetParticle.GetSide()!=-1)
    {
      backwardEnergy += targetMass;
      backwardCount --;
      forwardEnergy -= targetMass;
      forwardCount ++;
    }
  }

  for (i=0; i<vecLen; ++i) {
    if( vec[i]->GetSide() == -1 )
    {
      ++backwardCount;
      backwardEnergy -= vec[i]->GetMass()/GeV;
    } else {
      ++forwardCount;
      forwardEnergy -= vec[i]->GetMass()/GeV;
    }
  }

  // Check that sum of forward particle masses does not exceed forwardEnergy, 
  // and similarly for backward particles.  If so, move particles from one
  // hemisphere to another.

  if (backwardEnergy < 0.0) {
    for (i = 0; i < vecLen; ++i) {
      if (vec[i]->GetSide() == -1) { 
        backwardEnergy += vec[i]->GetMass()/GeV;
        --backwardCount;
        vec[i]->SetSide(1);
        forwardEnergy -= vec[i]->GetMass()/GeV;
        ++forwardCount;
        if (backwardEnergy > 0.0) break;
      }
    }
  }

  if (forwardEnergy < 0.0) {
    for (i = 0; i < vecLen; ++i) {
      if (vec[i]->GetSide() == 1) { 
        forwardEnergy += vec[i]->GetMass()/GeV;
        --forwardCount;
        vec[i]->SetSide(-1);
        backwardEnergy -= vec[i]->GetMass()/GeV;
        ++backwardCount;
        if (forwardEnergy > 0.0) break;
      }
    }
  }
 
  // Special cases for reactions near threshold

  // 1. There is only one secondary 
  if (forwardEnergy > 0.0 && backwardEnergy < 0.0) {
    forwardEnergy += backwardEnergy;
    backwardEnergy = 0;
  }

  // 2. Nuclear binding energy is large 
  if (forwardEnergy + backwardEnergy < 0.0) return false;


  // forwardEnergy now represents the total energy in the forward reaction 
  // hemisphere which is available for kinetic energy and particle creation.
  // Similarly for backwardEnergy.

  // Add particles from the intranuclear cascade.
  // nuclearExcitationCount = number of new secondaries produced by nuclear 
  // excitation
  // extraCount = number of nucleons within these new secondaries
  //
  // Note: eventually have to make sure that enough nucleons are available 
  // in the case of small target nuclei

  G4double xtarg;
  G4double a13 = G4Pow::GetInstance()->A13(atomicWeight);  // A**(1/3)
  if (centerofmassEnergy < (2.0+G4UniformRand()) )
    xtarg = afc * (a13-1.0) * (2.0*backwardCount+vecLen+2)/2.0;
  else
    xtarg = afc * (a13-1.0) * (2.0*backwardCount);

  if( xtarg <= 0.0 )xtarg = 0.01;
  G4int nuclearExcitationCount = G4Poisson( xtarg );
  // To do: try reducing nuclearExcitationCount with increasing energy 
  //        to simulate cut-off of cascade
  if(atomicWeight<1.0001) nuclearExcitationCount = 0;
  G4int extraNucleonCount = 0;
  G4double extraNucleonMass = 0.0;

  if (nuclearExcitationCount > 0) {
    const G4double nucsup[] = { 1.00, 0.7, 0.5, 0.4, 0.35, 0.3 };
    const G4double psup[] = { 3., 6., 20., 50., 100., 1000. };
    G4int momentumBin = 0;

    G4int loop = 0;
    G4ExceptionDescription ed;
    ed << " While count exceeded " << G4endl;
    while( (momentumBin < 6) &&
           (modifiedOriginal.GetTotalMomentum()/GeV > psup[momentumBin]) ) {  /* Loop checking, 01.09.2015, D.Wright */
      ++momentumBin;
      loop++;
      if (loop > 1000) {
        G4Exception("G4RPGFragmentation::ReactionStage()", "HAD_RPG_100", JustWarning, ed);
        break;
      }
    }

    momentumBin = std::min( 5, momentumBin );

    //  NOTE: in GENXPT, these new particles were given negative codes
    //        here I use  NewlyAdded = true  instead
    //

    for (i = 0; i < nuclearExcitationCount; ++i) {
      G4ReactionProduct * pVec = new G4ReactionProduct();
      if (G4UniformRand() < nucsup[momentumBin]) {

        if (G4UniformRand() > 1.0-atomicNumber/atomicWeight)
          pVec->SetDefinition( aProton );
        else
          pVec->SetDefinition( aNeutron );

        // nucleon comes from nucleus - 
	// do not subtract its mass from backward energy
        pVec->SetSide( -2 );                // -2 means backside nucleon
        ++extraNucleonCount;
        extraNucleonMass += pVec->GetMass()/GeV;
	// To do: Remove chosen nucleon from target nucleus
        pVec->SetNewlyAdded( true );
        vec.SetElement( vecLen++, pVec );    
        ++backwardCount;

      } else {

        G4double ran = G4UniformRand();
        if( ran < 0.3181 )
          pVec->SetDefinition( aPiPlus );
        else if( ran < 0.6819 )
          pVec->SetDefinition( aPiZero );
        else
          pVec->SetDefinition( aPiMinus );

        if (backwardEnergy > pVec->GetMass()/GeV) {
          backwardEnergy -= pVec->GetMass()/GeV;    // pion mass comes from free energy
          ++backwardCount;
          pVec->SetSide( -1 );         // backside particle, but not a nucleon
          pVec->SetNewlyAdded( true );
          vec.SetElement( vecLen++, pVec );
	}

	// To do: Change proton to neutron (or vice versa) in target nucleus depending 
	//        on pion charge
      }
    }
  }

  // Define initial state vectors for Lorentz transformations
  // The pseudoParticles have non-standard masses, hence the "pseudo"

  G4ReactionProduct pseudoParticle[8];
  for (i = 0; i < 8; ++i) pseudoParticle[i].SetZero();
    
  pseudoParticle[0].SetMass( mOriginal*GeV );
  pseudoParticle[0].SetMomentum( 0.0, 0.0, pOriginal*GeV );
  pseudoParticle[0].SetTotalEnergy(
   std::sqrt( pOriginal*pOriginal + mOriginal*mOriginal )*GeV );

  pseudoParticle[1].SetMass(protonMass);     // this could be targetMass
  pseudoParticle[1].SetTotalEnergy(protonMass);
    
  pseudoParticle[3].SetMass(protonMass*(1+extraNucleonCount) );
  pseudoParticle[3].SetTotalEnergy(protonMass*(1+extraNucleonCount) );
    
  pseudoParticle[2] = pseudoParticle[0] + pseudoParticle[1];
  pseudoParticle[3] = pseudoParticle[3] + pseudoParticle[0];
    
  pseudoParticle[0].Lorentz( pseudoParticle[0], pseudoParticle[2] );
  pseudoParticle[1].Lorentz( pseudoParticle[1], pseudoParticle[2] );
    
  // Main loop for 4-momentum generation
  // See Pitha-report (Aachen) for a detailed description of the method

  G4double aspar, pt, et, x, pp, pp1, wgt;
  G4int    innerCounter, outerCounter;
  G4bool   eliminateThisParticle, resetEnergies, constantCrossSection;
    
  G4double forwardKinetic = 0.0;
  G4double backwardKinetic = 0.0;

  // Process the secondary particles in reverse order
  // The incident particle is done after the secondaries
  // Nucleons, including the target, in the backward hemisphere are also 
  // done later

  G4int backwardNucleonCount = 0;  // number of nucleons in backward hemisphere
  G4double totalEnergy, kineticEnergy, vecMass;
  G4double phi;

  for (i = vecLen-1; i >= 0; --i) {

    if (vec[i]->GetNewlyAdded()) {         // added from intranuclear cascade
      if (vec[i]->GetSide() == -2) {       // its a nucleon
        if (backwardNucleonCount < 18) {
          if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
            for (G4int j = 0; j < vecLen; j++) delete vec[j];
            vecLen = 0;
            throw G4HadReentrentException(__FILE__, __LINE__,
            "G4RPGFragmentation::ReactionStage : a pion has been counted as a backward nucleon");
          }
          vec[i]->SetSide(-3);
          ++backwardNucleonCount;
          continue;     // Don't generate momenta for the first 17 backward 
                        // cascade nucleons. This gets done by the cluster 
	                // model later on.
        }
      }
    }

    // Set pt and phi values, they are changed somewhat in the iteration loop
    // Set mass parameter for lambda fragmentation model

    vecMass = vec[i]->GetMass()/GeV;
    G4double ran = -G4Log(1.0-G4UniformRand())/3.5;

    if (vec[i]->GetSide() == -2) {  // backward nucleon
      aspar = 0.20;
      pt = std::sqrt( std::pow( ran, 1.2 ) );

    } else {                        // not a backward nucleon
      if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
        aspar = 0.75;
        pt = std::sqrt( std::pow( ran, 1.7 ) );
      } else if (vec[i]->GetDefinition()->GetParticleSubType() == "kaon") {
        aspar = 0.70;
        pt = std::sqrt( std::pow( ran, 1.7 ) );
      } else {                        // vec[i] must be a baryon or ion 
        aspar = 0.65;
        pt = std::sqrt( std::pow( ran, 1.5 ) );
      }
    }

    pt = std::max( 0.001, pt );
    phi = G4UniformRand()*twopi;
    vec[i]->SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );
    if (vec[i]->GetSide() > 0)
      et = pseudoParticle[0].GetTotalEnergy()/GeV;
    else
      et = pseudoParticle[1].GetTotalEnergy()/GeV;

    //
    // Start of outer iteration loop
    //
    outerCounter = 0;
    eliminateThisParticle = true;
    resetEnergies = true;
    dndl[0] = 0.0;

    while (++outerCounter < 3) {  /* Loop checking, 01.09.2015, D.Wright */

      FragmentationIntegral(pt, et, aspar, vecMass);
      innerCounter = 0;
      vec[i]->SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );

      // Start of inner iteration loop

      while (++innerCounter < 7) {  /* Loop checking, 01.09.2015, D.Wright */

        ran = G4UniformRand()*dndl[19];
        l = 1;

        G4int loop = 0;
        G4ExceptionDescription ed;
        ed << " While count exceeded " << G4endl; 
        while( ( ran > dndl[l] ) && ( l < 19 ) ) { /* Loop checking, 01.09.2015, D.Wright */
          l++;
          loop++;
          if (loop > 1000) {
            G4Exception("G4RPGFragmentation::ReactionStage()", "HAD_RPG_100", JustWarning, ed);
            break;
          }
        }

        x = (G4double(l-1) + G4UniformRand())/19.;
        if (vec[i]->GetSide() < 0) x *= -1.;
        vec[i]->SetMomentum( x*et*GeV );              // set the z-momentum
        totalEnergy = std::sqrt( x*et*x*et + pt*pt + vecMass*vecMass );
        vec[i]->SetTotalEnergy( totalEnergy*GeV );
        kineticEnergy = vec[i]->GetKineticEnergy()/GeV;

        if (vec[i]->GetSide() > 0) {                  // forward side
          if( (forwardKinetic+kineticEnergy) < 0.95*forwardEnergy ) {
	    // Leave at least 5% of the forward free energy for the projectile

            pseudoParticle[4] = pseudoParticle[4] + (*vec[i]);
            forwardKinetic += kineticEnergy;
            outerCounter = 2;                  // leave outer loop
            eliminateThisParticle = false;     // don't eliminate this particle
            resetEnergies = false;
            break;                             // leave inner loop
          }
          if( innerCounter > 5 )break;         // leave inner loop
          if( backwardEnergy >= vecMass )      // switch sides
          {
            vec[i]->SetSide(-1);
            forwardEnergy += vecMass;
            backwardEnergy -= vecMass;
            ++backwardCount;
          }
        } else {                               // backward side
	  //  if (extraNucleonCount > 19) x = 0.999;
	  //  G4double xxx = 0.95+0.05*extraNucleonCount/20.0;
          //  DHW: I think above lines were meant to be as follows:
          G4double xxx = 0.999;
          if (extraNucleonCount < 20) xxx = 0.95+0.05*extraNucleonCount/20.0;

          if ((backwardKinetic+kineticEnergy) < xxx*backwardEnergy) {
            pseudoParticle[5] = pseudoParticle[5] + (*vec[i]);
            backwardKinetic += kineticEnergy;
            outerCounter = 2;                  // leave outer loop
            eliminateThisParticle = false;     // don't eliminate this particle
            resetEnergies = false;
            break;                             // leave inner loop
          }
          if (innerCounter > 5) break;         // leave inner loop
          if (forwardEnergy >= vecMass) {      // switch sides
            vec[i]->SetSide(1);
            forwardEnergy -= vecMass;
            backwardEnergy += vecMass;
            backwardCount--;
          }
        }
        G4ThreeVector momentum = vec[i]->GetMomentum();
        vec[i]->SetMomentum( momentum.x() * 0.9, momentum.y() * 0.9 );
        pt *= 0.9;
        dndl[19] *= 0.9;
      }                      // closes inner loop

      // If we get here, the inner loop has been done 6 times.
      // If necessary, reduce energies of the previously done particles if 
      // they are lighter than protons or are in the forward hemisphere.
      // Then continue with outer loop.

      if (resetEnergies)
        ReduceEnergiesOfSecondaries(i+1, forwardKinetic, backwardKinetic,
                                    vec, vecLen,
                                    pseudoParticle[4], pseudoParticle[5],
                                    pt);

    } // closes outer loop
              
    if (eliminateThisParticle && vec[i]->GetMayBeKilled()) {
      // not enough energy, eliminate this particle

      if (vec[i]->GetSide() > 0) {
        --forwardCount;
        forwardEnergy += vecMass;
      } else {
        --backwardCount;
        if (vec[i]->GetSide() == -2) {
          --extraNucleonCount;
          extraNucleonMass -= vecMass;
        } else {
          backwardEnergy += vecMass;
	}
      }

      for( G4int j=i; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
      G4ReactionProduct* temp = vec[vecLen-1];
      delete temp;
      // To do: modify target nucleus according to particle eliminated  

      if( --vecLen == 0 ){
        G4cout << " FALSE RETURN DUE TO ENERGY BALANCE " << G4endl;
        return false;
      }  // all the secondaries have been eliminated
    }
  } // closes main loop over secondaries

  // Re-balance forward and backward energy if possible and if necessary

  G4double forwardKEDiff = forwardEnergy - forwardKinetic;
  G4double backwardKEDiff = backwardEnergy - backwardKinetic;

  if (forwardKEDiff < 0.0 || backwardKEDiff < 0.0) {
    ReduceEnergiesOfSecondaries(0, forwardKinetic, backwardKinetic,
                                vec, vecLen,
                                pseudoParticle[4], pseudoParticle[5],
                                pt);
  
    forwardKEDiff = forwardEnergy - forwardKinetic;
    backwardKEDiff = backwardEnergy - backwardKinetic;
    if (backwardKEDiff < 0.0) {
      if (forwardKEDiff + backwardKEDiff > 0.0) {
        backwardEnergy = backwardKinetic;
        forwardEnergy += backwardKEDiff;
        forwardKEDiff = forwardEnergy - forwardKinetic;
        backwardKEDiff = 0.0;
      } else {
        G4cout << " False return due to insufficient backward energy " << G4endl; 
        return false;
      }
    }
      
    if (forwardKEDiff < 0.0) {
      if (forwardKEDiff + backwardKEDiff > 0.0) {
        forwardEnergy = forwardKinetic;
        backwardEnergy += forwardKEDiff;
        backwardKEDiff = backwardEnergy - backwardKinetic;
        forwardKEDiff = 0.0;
      } else {
        G4cout << " False return due to insufficient forward energy " << G4endl; 
        return false;
      }
    }
  }

  // Generate momentum for incident (current) particle, which was placed 
  // in the forward hemisphere. 
  // Set mass parameter for lambda fragmentation model.
  // Set pt and phi values, which are changed somewhat in the iteration loop

  G4double ran = -G4Log(1.0-G4UniformRand());
  if (currentParticle.GetDefinition()->GetParticleSubType() == "pi") {
    aspar = 0.60;
    pt = std::sqrt( std::pow( ran/6.0, 1.7 ) );
  } else if (currentParticle.GetDefinition()->GetParticleSubType() == "kaon") {
    aspar = 0.50;
    pt = std::sqrt( std::pow( ran/5.0, 1.4 ) );
  } else {
    aspar = 0.40;
    pt = std::sqrt( std::pow( ran/4.0, 1.2 ) );
  }

  phi = G4UniformRand()*twopi;
  currentParticle.SetMomentum(pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV);
  et = pseudoParticle[0].GetTotalEnergy()/GeV;
  dndl[0] = 0.0;
  vecMass = currentParticle.GetMass()/GeV;

  FragmentationIntegral(pt, et, aspar, vecMass);

  ran = G4UniformRand()*dndl[19];
  l = 1;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
  while( ( ran > dndl[l] ) && ( l < 19 ) ) { /* Loop checking, 01.09.2015, D.Wright */
    l++;
    loop++;
    if (loop > 1000) {
      G4Exception("G4RPGFragmentation::ReactionStage()", "HAD_RPG_100", JustWarning, ed);
      break;
    }
  }

  x = (G4double(l-1) + G4UniformRand())/19.;
  currentParticle.SetMomentum( x*et*GeV );             // set the z-momentum

  if (forwardEnergy < forwardKinetic) {
    totalEnergy = vecMass + 0.04*std::fabs(normal());
    G4cout << " Not enough forward energy: forwardEnergy = " 
           << forwardEnergy << " forwardKinetic = "  
           << forwardKinetic << " total energy left = " 
           << backwardKEDiff + forwardKEDiff << G4endl;
  } else {
    totalEnergy = vecMass + forwardEnergy - forwardKinetic;
    forwardKinetic = forwardEnergy;
  }
  currentParticle.SetTotalEnergy( totalEnergy*GeV );
  pp = std::sqrt(std::abs( totalEnergy*totalEnergy - vecMass*vecMass) )*GeV;
  pp1 = currentParticle.GetMomentum().mag();

  if (pp1 < 1.0e-6*GeV) {
    G4ThreeVector iso = Isotropic(pp);
    currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
  } else {
    currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
  }
  pseudoParticle[4] = pseudoParticle[4] + currentParticle;

  // Current particle now finished

  // Begin target particle

  if (backwardNucleonCount < 18) {
    targetParticle.SetSide(-3);
    ++backwardNucleonCount;

  } else {
    // Set pt and phi values, they are changed somewhat in the iteration loop
    // Set mass parameter for lambda fragmentation model

    vecMass = targetParticle.GetMass()/GeV;
    ran = -G4Log(1.0-G4UniformRand());
    aspar = 0.40;
    pt = std::max( 0.001, std::sqrt( std::pow( ran/4.0, 1.2 ) ) );
    phi = G4UniformRand()*twopi;
    targetParticle.SetMomentum(pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV);
    et = pseudoParticle[1].GetTotalEnergy()/GeV;
    outerCounter = 0;
    innerCounter = 0;
    G4bool marginalEnergy = true;
    dndl[0] = 0.0;
    G4double xxx = 0.999;
    if( extraNucleonCount < 20 ) xxx = 0.95+0.05*extraNucleonCount/20.0;
    G4ThreeVector momentum;

    while (++outerCounter < 4) {  /* Loop checking, 01.09.2015, D.Wright */
      FragmentationIntegral(pt, et, aspar, vecMass);

      for (innerCounter = 0; innerCounter < 6; innerCounter++) {
        ran = G4UniformRand()*dndl[19];
        l = 1;

        G4int loopa = 0;
        G4ExceptionDescription eda;
        eda << " While count exceeded " << G4endl; 
        while( ( ran > dndl[l] ) && ( l < 19 ) ) { /* Loop checking, 01.09.2015, D.Wright */
          l++;
          loopa++;
          if (loopa > 1000) {
            G4Exception("G4RPGFragmentation::ReactionStage()", "HAD_RPG_100", JustWarning, eda);
            break;
          }
        }

        x = -(G4double(l-1) + G4UniformRand())/19.;
        targetParticle.SetMomentum( x*et*GeV );        // set the z-momentum
        totalEnergy = std::sqrt(x*et*x*et + pt*pt + vecMass*vecMass);
        targetParticle.SetTotalEnergy( totalEnergy*GeV );

        if ((backwardKinetic+totalEnergy-vecMass) < xxx*backwardEnergy) {
          pseudoParticle[5] = pseudoParticle[5] + targetParticle;
          backwardKinetic += totalEnergy - vecMass;
          outerCounter = 3;                 // leave outer loop
          marginalEnergy = false;
          break;                            // leave inner loop
        }
        momentum = targetParticle.GetMomentum();
        targetParticle.SetMomentum(momentum.x() * 0.9, momentum.y() * 0.9);
        pt *= 0.9;
        dndl[19] *= 0.9;
      }
    }

    if (marginalEnergy) {
      G4cout << " Extra backward kinetic energy = "
             << 0.999*backwardEnergy - backwardKinetic << G4endl;
      totalEnergy = vecMass + 0.999*backwardEnergy - backwardKinetic;
      targetParticle.SetTotalEnergy(totalEnergy*GeV);
      pp = std::sqrt(std::abs(totalEnergy*totalEnergy - vecMass*vecMass) )*GeV;
      targetParticle.SetMomentum(momentum.x()/0.9, momentum.y()/0.9);
      pp1 = targetParticle.GetMomentum().mag();
      targetParticle.SetMomentum(targetParticle.GetMomentum() * pp/pp1 );
      pseudoParticle[5] = pseudoParticle[5] + targetParticle;
      backwardKinetic = 0.999*backwardEnergy;
    }

  }

  if (backwardEnergy < backwardKinetic) 
    G4cout << " Backward Edif = " << backwardEnergy - backwardKinetic << G4endl;
  if (forwardEnergy != forwardKinetic) 
    G4cout << " Forward Edif = " << forwardEnergy - forwardKinetic << G4endl;

  // Target particle finished.
  // Now produce backward nucleons with a cluster model
  // ps[2] = CM frame of projectile + target
  // ps[3] = sum of projectile + nucleon cluster in lab frame
  // ps[6] = proj + cluster 4-vector boosted into CM frame of proj + targ
  //         with secondaries, current and target particles subtracted
  //       = total 4-momentum of backward nucleon cluster  

  pseudoParticle[6].Lorentz( pseudoParticle[3], pseudoParticle[2] );
  pseudoParticle[6] = pseudoParticle[6] - pseudoParticle[4];
  pseudoParticle[6] = pseudoParticle[6] - pseudoParticle[5];

  if (backwardNucleonCount == 1) {   
    // Target particle is the only backward nucleon.  Give it the remainder
    // of the backward kinetic energy. 

    G4double ekin =
      std::min(backwardEnergy-backwardKinetic, centerofmassEnergy/2.0-protonMass/GeV);

    if( ekin < 0.04 )ekin = 0.04 * std::fabs( normal() );
    vecMass = targetParticle.GetMass()/GeV;
    totalEnergy = ekin + vecMass;
    targetParticle.SetTotalEnergy( totalEnergy*GeV );
    pp = std::sqrt(std::abs(totalEnergy*totalEnergy - vecMass*vecMass) )*GeV;
    pp1 = pseudoParticle[6].GetMomentum().mag();
    if (pp1 < 1.0e-6*GeV) { 
      G4ThreeVector iso = Isotropic(pp);
      targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
    } else {
      targetParticle.SetMomentum( pseudoParticle[6].GetMomentum() * (pp/pp1));
    }
    pseudoParticle[5] = pseudoParticle[5] + targetParticle;

  } else if (backwardNucleonCount > 1) {
    // Share remaining energy with up to 17 backward nucleons

    G4int tempCount = 5;
    if (backwardNucleonCount < 5) tempCount = backwardNucleonCount;
    tempCount -= 2;

    G4double clusterMass = 0.;
    if (targetParticle.GetSide() == -3) 
      clusterMass = targetParticle.GetMass()/GeV;
    for (i = 0; i < vecLen; ++i)
      if (vec[i]->GetSide() == -3) clusterMass += vec[i]->GetMass()/GeV;
    clusterMass += backwardEnergy - backwardKinetic;

    totalEnergy = pseudoParticle[6].GetTotalEnergy()/GeV;
    pseudoParticle[6].SetMass(clusterMass*GeV);

    pp = std::sqrt(std::abs(totalEnergy*totalEnergy - 
                            clusterMass*clusterMass) )*GeV;
    pp1 = pseudoParticle[6].GetMomentum().mag();
    if (pp1 < 1.0e-6*GeV) {
      G4ThreeVector iso = Isotropic(pp);
      pseudoParticle[6].SetMomentum(iso.x(), iso.y(), iso.z());
    } else {
      pseudoParticle[6].SetMomentum(pseudoParticle[6].GetMomentum() * (-pp/pp1));
    }

    std::vector<G4ReactionProduct*> tempList;  // ptrs to backward nucleons
    if (targetParticle.GetSide() == -3) tempList.push_back(&targetParticle);
    for (i = 0; i < vecLen; ++i)
      if (vec[i]->GetSide() == -3) tempList.push_back(vec[i]);

    constantCrossSection = true;

    if (tempList.size() > 1) {
      G4int n_entry = 0;
      wgt = GenerateNBodyEventT(pseudoParticle[6].GetMass(), 
                                constantCrossSection, tempList);

      if (targetParticle.GetSide() == -3) {
        targetParticle = *tempList[0];
        targetParticle.Lorentz(targetParticle, pseudoParticle[6]);
        n_entry++;
      }

      for (i = 0; i < vecLen; ++i) {
        if (vec[i]->GetSide() == -3) {
          *vec[i] = *tempList[n_entry];
          vec[i]->Lorentz(*vec[i], pseudoParticle[6]);
          n_entry++;
        }
      }
    }
  } else return false;

  if (vecLen == 0) return false;  // all the secondaries have been eliminated

  // Lorentz transformation to lab frame

  currentParticle.Lorentz( currentParticle, pseudoParticle[1] );
  targetParticle.Lorentz( targetParticle, pseudoParticle[1] );    
  for (i = 0; i < vecLen; ++i) vec[i]->Lorentz(*vec[i], pseudoParticle[1]);

  // Set leading strange particle flag.
  // leadFlag will be true if original particle and incident particle are
  // both strange, in which case the incident particle becomes the leading 
  // particle. 
  // leadFlag will also be true if the target particle is strange, but the
  // incident particle is not, in which case the target particle becomes the 
  // leading particle.

  G4bool leadingStrangeParticleHasChanged = true;
  if (leadFlag)
  {
    if (currentParticle.GetDefinition() == leadingStrangeParticle.GetDefinition())
      leadingStrangeParticleHasChanged = false;
    if (leadingStrangeParticleHasChanged &&
        (targetParticle.GetDefinition() == leadingStrangeParticle.GetDefinition()) )
      leadingStrangeParticleHasChanged = false;
    if( leadingStrangeParticleHasChanged )
    {
      for( i=0; i<vecLen; i++ )
      {
        if( vec[i]->GetDefinition() == leadingStrangeParticle.GetDefinition() )
        {
          leadingStrangeParticleHasChanged = false;
          break;
        }
      }
    }
    if( leadingStrangeParticleHasChanged )
    {
      G4bool leadTest = 
        (leadingStrangeParticle.GetDefinition()->GetParticleSubType() == "kaon" ||
         leadingStrangeParticle.GetDefinition()->GetParticleSubType() == "pi");
      G4bool targetTest =
        (targetParticle.GetDefinition()->GetParticleSubType() == "kaon" ||
         targetParticle.GetDefinition()->GetParticleSubType() == "pi");
        
      // following modified by JLC 22-Oct-97
          
      if( (leadTest&&targetTest) || !(leadTest||targetTest) ) // both true or both false
      {
        targetParticle.SetDefinitionAndUpdateE( leadingStrangeParticle.GetDefinition() );
        targetHasChanged = true;
      }
      else
      {
        currentParticle.SetDefinitionAndUpdateE( leadingStrangeParticle.GetDefinition() );
        incidentHasChanged = false;
      }
    }
  }   // end of if( leadFlag )

  // Get number of final state nucleons and nucleons remaining in
  // target nucleus
    
  std::pair<G4int, G4int> finalStateNucleons = 
    GetFinalStateNucleons(originalTarget, vec, vecLen);

  G4int protonsInFinalState = finalStateNucleons.first;
  G4int neutronsInFinalState = finalStateNucleons.second;

  G4int numberofFinalStateNucleons = 
    protonsInFinalState + neutronsInFinalState;

  if (currentParticle.GetDefinition()->GetBaryonNumber() == 1 &&
      targetParticle.GetDefinition()->GetBaryonNumber() == 1 &&
      originalIncident->GetDefinition()->GetPDGMass() < 
                                 G4Lambda::Lambda()->GetPDGMass())
    numberofFinalStateNucleons++;

  numberofFinalStateNucleons = std::max(1, numberofFinalStateNucleons);

  G4int PinNucleus = std::max(0, 
    G4int(targetNucleus.GetZ_asInt()) - protonsInFinalState);
  G4int NinNucleus = std::max(0,
    G4int(targetNucleus.GetA_asInt()-targetNucleus.GetZ_asInt()) - neutronsInFinalState);

  pseudoParticle[3].SetMomentum( 0.0, 0.0, pOriginal*GeV );
  pseudoParticle[3].SetMass( mOriginal*GeV );
  pseudoParticle[3].SetTotalEnergy(
   std::sqrt( pOriginal*pOriginal + mOriginal*mOriginal )*GeV );
    
  const G4ParticleDefinition* aOrgDef = modifiedOriginal.GetDefinition();
  G4int diff = 0;
  if(aOrgDef == G4Proton::Proton() || aOrgDef == G4Neutron::Neutron() )  diff = 1;
  if(numberofFinalStateNucleons == 1) diff = 0;
  pseudoParticle[4].SetMomentum( 0.0, 0.0, 0.0 );
  pseudoParticle[4].SetMass( protonMass*(numberofFinalStateNucleons-diff) );
  pseudoParticle[4].SetTotalEnergy( protonMass*(numberofFinalStateNucleons-diff) );
    
  G4double theoreticalKinetic = 
    pseudoParticle[3].GetTotalEnergy() + pseudoParticle[4].GetTotalEnergy() -
    currentParticle.GetMass() - targetParticle.GetMass();
  for (i = 0; i < vecLen; ++i) theoreticalKinetic -= vec[i]->GetMass();
    
  G4double simulatedKinetic = 
    currentParticle.GetKineticEnergy() + targetParticle.GetKineticEnergy();
  for (i = 0; i < vecLen; ++i) 
    simulatedKinetic += vec[i]->GetKineticEnergy();
    
  pseudoParticle[5] = pseudoParticle[3] + pseudoParticle[4];
  pseudoParticle[3].Lorentz( pseudoParticle[3], pseudoParticle[5] );
  pseudoParticle[4].Lorentz( pseudoParticle[4], pseudoParticle[5] );
    
  pseudoParticle[7].SetZero();
  pseudoParticle[7] = pseudoParticle[7] + currentParticle;
  pseudoParticle[7] = pseudoParticle[7] + targetParticle;
  for (i = 0; i < vecLen; ++i)
    pseudoParticle[7] = pseudoParticle[7] + *vec[i];

    /*
    // This code does not appear to do anything to the energy or angle spectra
    if( vecLen <= 16 && vecLen > 0 )
    {
      // must create a new set of ReactionProducts here because GenerateNBody will
      // modify the momenta for the particles, and we don't want to do this
      //
      G4ReactionProduct tempR[130];
      tempR[0] = currentParticle;
      tempR[1] = targetParticle;
      for( i=0; i<vecLen; ++i )tempR[i+2] = *vec[i];
      G4FastVector<G4ReactionProduct,256> tempV1;
      tempV1.Initialize( vecLen+2 );
      G4int tempLen1 = 0;
      for( i=0; i<vecLen+2; ++i )tempV1.SetElement( tempLen1++, &tempR[i] );
      constantCrossSection = true;

      wgt = GenerateNBodyEvent(pseudoParticle[3].GetTotalEnergy() +
                               pseudoParticle[4].GetTotalEnergy(),
                               constantCrossSection, tempV1, tempLen1);
      if (wgt == -1) {
        G4double Qvalue = 0;
        for (i = 0; i < tempLen1; i++) Qvalue += tempV1[i]->GetMass();
        wgt = GenerateNBodyEvent(Qvalue,
                                 constantCrossSection, tempV1, tempLen1);
      }
      if(wgt>-.5)
      {
        theoreticalKinetic = 0.0;
        for( i=0; i<tempLen1; ++i )
        {
          pseudoParticle[6].Lorentz( *tempV1[i], pseudoParticle[4] );
          theoreticalKinetic += pseudoParticle[6].GetKineticEnergy();
        }
      }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    }
    */

  //
  // Make sure that the kinetic energies are correct
  //

  if (simulatedKinetic != 0.0) {
    wgt = theoreticalKinetic/simulatedKinetic;
    theoreticalKinetic = currentParticle.GetKineticEnergy() * wgt;
    simulatedKinetic = theoreticalKinetic;
    currentParticle.SetKineticEnergy(theoreticalKinetic);
    pp = currentParticle.GetTotalMomentum();
    pp1 = currentParticle.GetMomentum().mag();
    if (pp1 < 1.0e-6*GeV) {
      G4ThreeVector iso = Isotropic(pp);
      currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
    } else {
      currentParticle.SetMomentum(currentParticle.GetMomentum() * (pp/pp1));
    }

    theoreticalKinetic = targetParticle.GetKineticEnergy() * wgt;
    targetParticle.SetKineticEnergy(theoreticalKinetic);
    simulatedKinetic += theoreticalKinetic;
    pp = targetParticle.GetTotalMomentum();
    pp1 = targetParticle.GetMomentum().mag();

    if (pp1 < 1.0e-6*GeV) {
      G4ThreeVector iso = Isotropic(pp);
      targetParticle.SetMomentum(iso.x(), iso.y(), iso.z() );
    } else {
      targetParticle.SetMomentum(targetParticle.GetMomentum() * (pp/pp1) );
    }

    for (i = 0; i < vecLen; ++i ) {
      theoreticalKinetic = vec[i]->GetKineticEnergy() * wgt;
      simulatedKinetic += theoreticalKinetic;
      vec[i]->SetKineticEnergy(theoreticalKinetic);
      pp = vec[i]->GetTotalMomentum();
      pp1 = vec[i]->GetMomentum().mag();
      if( pp1 < 1.0e-6*GeV ) {
        G4ThreeVector iso = Isotropic(pp);
        vec[i]->SetMomentum(iso.x(), iso.y(), iso.z() );
      } else {
        vec[i]->SetMomentum(vec[i]->GetMomentum() * (pp/pp1) );
      }
    }
  }

  //    Rotate(numberofFinalStateNucleons, pseudoParticle[3].GetMomentum(),
  //           modifiedOriginal, originalIncident, targetNucleus,
  //           currentParticle, targetParticle, vec, vecLen );

  // Add black track particles
  // the total number of particles produced is restricted to 198
  // this may have influence on very high energies

  if( atomicWeight >= 1.5 )
  {
    // npnb is number of proton/neutron black track particles
    // ndta is the number of deuterons, tritons, and alphas produced
    // epnb is the kinetic energy available for proton/neutron black track 
    //    particles
    // edta is the kinetic energy available for deuteron/triton/alpha particles

    G4int npnb = 0;
    G4int ndta = 0;
      
    G4double epnb, edta;
    if (veryForward) {
      epnb = targetNucleus.GetAnnihilationPNBlackTrackEnergy();
      edta = targetNucleus.GetAnnihilationDTABlackTrackEnergy();
    } else {
      epnb = targetNucleus.GetPNBlackTrackEnergy();
      edta = targetNucleus.GetDTABlackTrackEnergy();
    }
      /*
      G4ReactionProduct* fudge = new G4ReactionProduct();
      fudge->SetDefinition( aProton );
      G4double TT = epnb + edta;
      G4double MM = fudge->GetMass()/GeV;
      fudge->SetTotalEnergy(MM*GeV + TT*GeV);
      G4double pzz = std::sqrt(TT*(TT + 2.*MM));
      fudge->SetMomentum(0.0, 0.0, pzz*GeV);
      vec.SetElement(vecLen++, fudge);    
      // G4cout << " Fudge = " << vec[vecLen-1]->GetKineticEnergy()/GeV 
                << G4endl;
      */

    const G4double pnCutOff = 0.001;
    const G4double dtaCutOff = 0.001;
    //    const G4double kineticMinimum = 1.e-6;
    //    const G4double kineticFactor = -0.010;
    //    G4double sprob = 0.0;  // sprob = probability of self-absorption in 
                           // heavy molecules
    //  Not currently used (DHW 9 June 2008)  const G4double ekIncident = originalIncident->GetKineticEnergy()/GeV;
    //    if (ekIncident >= 5.0) sprob = std::min(1.0, 0.6*std::log(ekIncident-4.0));
    if (epnb > pnCutOff)
    {
      npnb = G4Poisson((1.5+1.25*numberofFinalStateNucleons)*epnb/(epnb+edta));
      if (numberofFinalStateNucleons + npnb > atomicWeight)
        npnb = G4int(atomicWeight+0.00001 - numberofFinalStateNucleons);
      npnb = std::min( npnb, 127-vecLen );
    }
    if( edta >= dtaCutOff )
    {
      ndta = G4Poisson((1.5+1.25*numberofFinalStateNucleons)*edta/(epnb+edta));
      ndta = std::min( ndta, 127-vecLen );
    }
    if (npnb == 0 && ndta == 0) npnb = 1;

    AddBlackTrackParticles(epnb, npnb, edta, ndta, modifiedOriginal,
                           PinNucleus, NinNucleus, targetNucleus,
                           vec, vecLen);
  }

  //  if( centerofmassEnergy <= (4.0+G4UniformRand()) )
  //    MomentumCheck( modifiedOriginal, currentParticle, targetParticle, 
  //                     vec, vecLen );
  //
  //  calculate time delay for nuclear reactions
  //

  if( (atomicWeight >= 1.5) && (atomicWeight <= 230.0) && (ekOriginal <= 0.2) )
    currentParticle.SetTOF( 
         1.0-500.0*G4Exp(-ekOriginal/0.04)*G4Log(G4UniformRand()) );
  else
    currentParticle.SetTOF( 1.0 );
  return true;

}


void G4RPGFragmentation::
ReduceEnergiesOfSecondaries(G4int startingIndex,
                            G4double& forwardKinetic,
                            G4double& backwardKinetic,
                            G4FastVector<G4ReactionProduct,256>& vec,
                            G4int& vecLen,
                            G4ReactionProduct& forwardPseudoParticle,
                            G4ReactionProduct& backwardPseudoParticle,
                            G4double& pt)
{
  // Reduce energies and pt of secondaries in order to maintain 
  // energy conservation

  G4double totalEnergy;
  G4double pp;
  G4double pp1;
  G4double px;
  G4double py;
  G4double mass;
  G4ReactionProduct* pVec; 
  G4int i;

  forwardKinetic = 0.0;
  backwardKinetic = 0.0;
  forwardPseudoParticle.SetZero();
  backwardPseudoParticle.SetZero();

  for (i = startingIndex; i < vecLen; i++) {
    pVec = vec[i];
    if (pVec->GetSide() != -3) {
      mass = pVec->GetMass();
      totalEnergy = 0.95*pVec->GetTotalEnergy() + 0.05*mass;
      pVec->SetTotalEnergy(totalEnergy);
      pp = std::sqrt( std::abs( totalEnergy*totalEnergy - mass*mass ) );
      pp1 = pVec->GetMomentum().mag();
      if (pp1 < 1.0e-6*GeV) {
        G4ThreeVector iso = Isotropic(pp);
        pVec->SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        pVec->SetMomentum(pVec->GetMomentum() * (pp/pp1) );
      }

      px = pVec->GetMomentum().x();
      py = pVec->GetMomentum().y();
      pt = std::max(1.0, std::sqrt( px*px + py*py ) )/GeV;
      if (pVec->GetSide() > 0) {
        forwardKinetic += pVec->GetKineticEnergy()/GeV;
        forwardPseudoParticle = forwardPseudoParticle + (*pVec);
      } else {
        backwardKinetic += pVec->GetKineticEnergy()/GeV;
        backwardPseudoParticle = backwardPseudoParticle + (*pVec);
      }
    }
  }
}

 /* end of file */
