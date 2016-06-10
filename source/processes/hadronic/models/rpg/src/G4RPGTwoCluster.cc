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
// $Id: G4RPGTwoCluster.cc 94214 2015-11-09 08:18:05Z gcosmo $
//

#include <iostream>
#include <signal.h>

#include "G4RPGTwoCluster.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4HadReentrentException.hh"

G4RPGTwoCluster::G4RPGTwoCluster()
  : G4RPGReaction() {}


G4bool G4RPGTwoCluster::
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
  // Derived from H. Fesefeldt's FORTRAN code TWOCLU
  //
  // A simple two cluster model is used to generate x- and pt- values for 
  // incident, target, and all secondary particles. 
  // This should be sufficient for low energy interactions.

  G4int i;
  G4ParticleDefinition* aProton = G4Proton::Proton();
  G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
  G4ParticleDefinition* aPiPlus = G4PionPlus::PionPlus();
  G4ParticleDefinition* aPiMinus = G4PionMinus::PionMinus();
  G4ParticleDefinition* aPiZero = G4PionZero::PionZero();
  G4bool veryForward = false;

  const G4double protonMass = aProton->GetPDGMass()/MeV;
  const G4double ekOriginal = modifiedOriginal.GetKineticEnergy()/GeV;
  const G4double etOriginal = modifiedOriginal.GetTotalEnergy()/GeV;
  const G4double mOriginal = modifiedOriginal.GetMass()/GeV;
  const G4double pOriginal = modifiedOriginal.GetMomentum().mag()/GeV;
  G4double targetMass = targetParticle.GetDefinition()->GetPDGMass()/GeV;
  G4double centerofmassEnergy = std::sqrt(mOriginal*mOriginal +
                                          targetMass*targetMass +
                                          2.0*targetMass*etOriginal);  // GeV
  G4double currentMass = currentParticle.GetMass()/GeV;
  targetMass = targetParticle.GetMass()/GeV;

  if (currentMass == 0.0 && targetMass == 0.0) {
    G4double ek = currentParticle.GetKineticEnergy();
    G4ThreeVector mom = currentParticle.GetMomentum();
    currentParticle = *vec[0];
    targetParticle = *vec[1];
    for (i = 0; i < (vecLen-2); ++i) *vec[i] = *vec[i+2];
    if (vecLen < 2) {
      for (G4int j = 0; j < vecLen; j++) delete vec[j];
      vecLen = 0;
      throw G4HadReentrentException(__FILE__, __LINE__,
      "G4RPGTwoCluster::ReactionStage : Negative number of particles");
    }
    delete vec[vecLen-1];
    delete vec[vecLen-2];
    vecLen -= 2;
    currentMass = currentParticle.GetMass()/GeV;
    targetMass = targetParticle.GetMass()/GeV;
    incidentHasChanged = true;
    targetHasChanged = true;
    currentParticle.SetKineticEnergy(ek);
    currentParticle.SetMomentum(mom);
    veryForward = true;
  }

  const G4double atomicWeight = targetNucleus.GetA_asInt();
  const G4double atomicNumber = targetNucleus.GetZ_asInt();

  // particles have been distributed in forward and backward hemispheres
  // in center of mass system of the hadron nucleon interaction

  // Incident particle always in forward hemisphere

  G4int forwardCount = 1;        // number of particles in forward hemisphere
  currentParticle.SetSide( 1 );
  G4double forwardMass = currentParticle.GetMass()/GeV;
  G4double cMass = forwardMass;
    
  // Target particle always in backward hemisphere
  G4int backwardCount = 1;       // number of particles in backward hemisphere
  targetParticle.SetSide( -1 );
  G4double backwardMass = targetParticle.GetMass()/GeV;
  G4double bMass = backwardMass;

  //  G4int backwardNucleonCount = 1;  // number of nucleons in backward hemisphere
  for (i = 0; i < vecLen; ++i) {
    if (vec[i]->GetSide() < 0) vec[i]->SetSide(-1);   // to take care of 
    // case where vec has been preprocessed by GenerateXandPt
    // and some of them have been set to -2 or -3
    if (vec[i]->GetSide() == -1) {
      ++backwardCount;
      backwardMass += vec[i]->GetMass()/GeV;
    } else {
      ++forwardCount;
      forwardMass += vec[i]->GetMass()/GeV;
    }
  }

  // Add nucleons and some pions from intra-nuclear cascade
  G4double term1 = G4Log(centerofmassEnergy*centerofmassEnergy);
  if(term1 < 0) term1 = 0.0001; // making sure xtarg<0;
  const G4double afc = 0.312 + 0.2 * G4Log(term1);
  G4double xtarg;
  G4double a13 = G4Pow::GetInstance()->A13(atomicWeight);  // A**(1/3)
  if( centerofmassEnergy < 2.0+G4UniformRand() )        // added +2 below, JLC 4Jul97
    xtarg = afc * (a13-1.0) * (2*backwardCount+vecLen+2)/2.0;
  else
    xtarg = afc * (a13-1.0) * (2*backwardCount);

  if( xtarg <= 0.0 )xtarg = 0.01;
  G4int nuclearExcitationCount = G4Poisson( xtarg );

  if(atomicWeight<1.0001) nuclearExcitationCount = 0;
  //  G4int extraNucleonCount = 0;
  //  G4double extraMass = 0.0;
  //  G4double extraNucleonMass = 0.0;
  if( nuclearExcitationCount > 0 )
  {
    G4int momentumBin = std::min( 4, G4int(pOriginal/3.0) );     
    const G4double nucsup[] = { 1.0, 0.8, 0.6, 0.5, 0.4 };
    //
    //  NOTE: in TWOCLU, these new particles were given negative codes
    //        here we use  NewlyAdded = true  instead
    //
    for( i=0; i<nuclearExcitationCount; ++i )
    {
      G4ReactionProduct* pVec = new G4ReactionProduct();
      if( G4UniformRand() < nucsup[momentumBin] )  // add proton or neutron
      {
        if( G4UniformRand() > 1.0-atomicNumber/atomicWeight )
          pVec->SetDefinition( aProton );
        else
          pVec->SetDefinition( aNeutron );
	// Not used        ++backwardNucleonCount;
	// Not used        ++extraNucleonCount;
	// Not used        extraNucleonMass += pVec->GetMass()/GeV;
      }
      else
      {                                       // add a pion
        G4double ran = G4UniformRand();
        if( ran < 0.3181 )
          pVec->SetDefinition( aPiPlus );
        else if( ran < 0.6819 )
          pVec->SetDefinition( aPiZero );
        else
          pVec->SetDefinition( aPiMinus );

	// DHW: add following two lines to correct energy balance
	//        ++backwardCount;
	//        backwardMass += pVec->GetMass()/GeV;
      }
      pVec->SetSide( -2 );    // backside particle
      // Not used     extraMass += pVec->GetMass()/GeV;
      pVec->SetNewlyAdded( true );
      vec.SetElement( vecLen++, pVec );
    }
  }

  // Masses of particles added from cascade not included in energy balance.
  // That's correct for nucleons from the intra-nuclear cascade but not for 
  // pions from the cascade.
 
  G4double forwardEnergy = (centerofmassEnergy-cMass-bMass)/2.0 +cMass - forwardMass;
  G4double backwardEnergy = (centerofmassEnergy-cMass-bMass)/2.0 +bMass - backwardMass;
  G4double eAvailable = centerofmassEnergy - (forwardMass+backwardMass);
  G4bool secondaryDeleted;
  G4double pMass;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
 // must eliminate a particle
  while( eAvailable <= 0.0 ) { /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 1000) {
      G4Exception("G4RPGTwoCluster::ReactionStage()", "HAD_RPG_100", JustWarning, ed);
      break;
    }

    secondaryDeleted = false;
    for( i=(vecLen-1); i>=0; --i )
    {
      if( vec[i]->GetSide() == 1 && vec[i]->GetMayBeKilled())
      {
        pMass = vec[i]->GetMass()/GeV;
        for( G4int j=i; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];     // shift up
        --forwardCount;
        forwardEnergy += pMass;
        forwardMass -= pMass;
        secondaryDeleted = true;
        break;
      }
      else if( vec[i]->GetSide() == -1 && vec[i]->GetMayBeKilled())
      {
        pMass = vec[i]->GetMass()/GeV;
        for( G4int j=i; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
        --backwardCount;
        backwardEnergy += pMass;
        backwardMass -= pMass;
        secondaryDeleted = true;
        break;
      }
    }  // breaks go down to here

    if( secondaryDeleted )
    {
      delete vec[vecLen-1];
      --vecLen;
        // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    }
    else
    {
      if( vecLen == 0 ) return false;  // all secondaries have been eliminated
      if( targetParticle.GetSide() == -1 )
      {
        pMass = targetParticle.GetMass()/GeV;
        targetParticle = *vec[0];
        for( G4int j=0; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
        --backwardCount;
        backwardEnergy += pMass;
        backwardMass -= pMass;
        secondaryDeleted = true;
      }
      else if( targetParticle.GetSide() == 1 )
      {
        pMass = targetParticle.GetMass()/GeV;
        targetParticle = *vec[0];
        for( G4int j=0; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
        --forwardCount;
        forwardEnergy += pMass;
        forwardMass -= pMass;
        secondaryDeleted = true;
      }

      if( secondaryDeleted )
      {
        delete vec[vecLen-1];
        --vecLen;
      }
      else
      {
        if( currentParticle.GetSide() == -1 )
        {
          pMass = currentParticle.GetMass()/GeV;
          currentParticle = *vec[0];
          for( G4int j=0; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
          --backwardCount;
          backwardEnergy += pMass;
          backwardMass -= pMass;
          secondaryDeleted = true;
        }
        else if( currentParticle.GetSide() == 1 )
        {
          pMass = currentParticle.GetMass()/GeV;
          currentParticle = *vec[0];
          for( G4int j=0; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
          --forwardCount;
          forwardEnergy += pMass;
          forwardMass -= pMass;
          secondaryDeleted = true;
        }
        if( secondaryDeleted )
        {
          delete vec[vecLen-1];
          --vecLen;
        }
        else break;

      }  // secondary not deleted 
    }  // secondary not deleted

    eAvailable = centerofmassEnergy - (forwardMass+backwardMass);
  } // while

  //
  // This is the start of the TwoCluster function
  // Choose multi-particle resonance masses by sampling 
  //    P(M) = gc[g(M-M0)]**(c-1) *exp[-(g(M-M0))**c] 
  // for M > M0
  //
  // Use for the forward and backward clusters, but not
  // the cascade cluster

  const G4double cpar[] = { 1.60, 1.35, 1.15, 1.10 };
  const G4double gpar[] = { 2.60, 1.80, 1.30, 1.20 };
  G4int ntc = 0;

  if (forwardCount < 1 || backwardCount < 1) return false;  // array bounds protection

  G4double rmc = forwardMass;
  if (forwardCount > 1) {
    ntc = std::min(3,forwardCount-2);
    rmc += std::pow(-G4Log(1.0-G4UniformRand()),1./cpar[ntc])/gpar[ntc];
  }
  G4double rmd = backwardMass;
  if( backwardCount > 1 ) {
    ntc = std::min(3,backwardCount-2);
    rmd += std::pow(-G4Log(1.0-G4UniformRand()),1./cpar[ntc])/gpar[ntc];
  }

  loop = 0;
  G4ExceptionDescription eda;
  eda << " While count exceeded " << G4endl;
  while( rmc+rmd > centerofmassEnergy ) { /* Loop checking, 01.09.2015, D.Wright */
    loop++;
    if (loop > 1000) {
      G4Exception("G4RPGTwoCluster::ReactionStage()", "HAD_RPG_100", JustWarning, eda);
      break;
    }
 
    if( (rmc <= forwardMass) && (rmd <= backwardMass) )
    {
      G4double temp = 0.999*centerofmassEnergy/(rmc+rmd);
      rmc *= temp;
      rmd *= temp;
    }
    else
    {
      rmc = 0.1*forwardMass + 0.9*rmc;
      rmd = 0.1*backwardMass + 0.9*rmd;
    }
  }

  G4ReactionProduct pseudoParticle[8];
  for( i=0; i<8; ++i )pseudoParticle[i].SetZero();
    
  pseudoParticle[1].SetMass( mOriginal*GeV );
  pseudoParticle[1].SetTotalEnergy( etOriginal*GeV );
  pseudoParticle[1].SetMomentum( 0.0, 0.0, pOriginal*GeV );
    
  pseudoParticle[2].SetMass( protonMass*MeV );
  pseudoParticle[2].SetTotalEnergy( protonMass*MeV );
  pseudoParticle[2].SetMomentum( 0.0, 0.0, 0.0 );
  //
  //  transform into center of mass system
  //
  pseudoParticle[0] = pseudoParticle[1] + pseudoParticle[2];
  pseudoParticle[1].Lorentz( pseudoParticle[1], pseudoParticle[0] );
  pseudoParticle[2].Lorentz( pseudoParticle[2], pseudoParticle[0] );

  // Calculate cm momentum for forward and backward masses
  // W = sqrt(pf*pf + rmc*rmc) + sqrt(pf*pf + rmd*rmd)
  // Solve for pf

  const G4double pfMin = 0.0001;
  G4double pf = (centerofmassEnergy*centerofmassEnergy+rmd*rmd-rmc*rmc);
  pf *= pf;
  pf -= 4*centerofmassEnergy*centerofmassEnergy*rmd*rmd;
  pf = std::sqrt( std::max(pf,pfMin) )/(2.0*centerofmassEnergy);
  //
  //  set final state masses and energies in centre of mass system
  //
  pseudoParticle[3].SetMass( rmc*GeV );
  pseudoParticle[3].SetTotalEnergy( std::sqrt(pf*pf+rmc*rmc)*GeV );
    
  pseudoParticle[4].SetMass( rmd*GeV );
  pseudoParticle[4].SetTotalEnergy( std::sqrt(pf*pf+rmd*rmd)*GeV );

  //
  // Get cm scattering angle by sampling t from tmin to tmax
  //
  const G4double bMin = 0.01;
  const G4double b1 = 4.0;
  const G4double b2 = 1.6;
  G4double pin = pseudoParticle[1].GetMomentum().mag()/GeV;
  G4double dtb = 4.0*pin*pf*std::max( bMin, b1+b2*G4Log(pOriginal) );
  G4double factor = 1.0 - G4Exp(-dtb);
  G4double costheta = 1.0 + 2.0*G4Log(1.0 - G4UniformRand()*factor) / dtb;

  costheta = std::max(-1.0, std::min(1.0, costheta));
  G4double sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  G4double phi = G4UniformRand() * twopi;
  //
  // calculate final state momenta in centre of mass system
  //
  pseudoParticle[3].SetMomentum( pf*sintheta*std::cos(phi)*GeV,
                                 pf*sintheta*std::sin(phi)*GeV,
                                 pf*costheta*GeV );
  pseudoParticle[4].SetMomentum( -pseudoParticle[3].GetMomentum());

  // Backward cluster of nucleons and pions from intra-nuclear cascade
  // Set up in lab system and transform to cms

  G4double pp, pp1;
  if( nuclearExcitationCount > 0 )
  {
    const G4double ga = 1.2;
    G4double ekit1 = 0.04;
    G4double ekit2 = 0.6;   // Max KE of cascade particle
    if( ekOriginal <= 5.0 )
    {
      ekit1 *= ekOriginal*ekOriginal/25.0;
      ekit2 *= ekOriginal*ekOriginal/25.0;
    }
    G4double scale = std::pow(ekit2/ekit1, 1.0-ga) - 1.0;
    for( i=0; i<vecLen; ++i )
    {
      if( vec[i]->GetSide() == -2 )
      {
        G4double kineticE = ekit1*std::pow((1.0 + G4UniformRand()*scale), 1.0/(1.0-ga) );
        vec[i]->SetKineticEnergy( kineticE*GeV );
        G4double vMass = vec[i]->GetMass()/MeV;
        G4double totalE = kineticE*GeV + vMass;
        pp = std::sqrt( std::abs(totalE*totalE-vMass*vMass) );
        G4double cost = std::min( 1.0, std::max( -1.0, G4Log(2.23*G4UniformRand()+0.383)/0.96 ) );
        G4double sint = std::sqrt(1.0-cost*cost);
        phi = twopi*G4UniformRand();
        vec[i]->SetMomentum( pp*sint*std::cos(phi)*MeV,
                             pp*sint*std::sin(phi)*MeV,
                             pp*cost*MeV );
        vec[i]->Lorentz( *vec[i], pseudoParticle[0] );
      }
    }
  }

  //
  // Fragmentation of forward and backward clusters
  //

  currentParticle.SetMomentum( pseudoParticle[3].GetMomentum() );
  currentParticle.SetTotalEnergy( pseudoParticle[3].GetTotalEnergy() );
    
  targetParticle.SetMomentum( pseudoParticle[4].GetMomentum() );
  targetParticle.SetTotalEnergy( pseudoParticle[4].GetTotalEnergy() );
    
  pseudoParticle[5].SetMomentum( pseudoParticle[3].GetMomentum() * (-1.0) );
  pseudoParticle[5].SetMass( pseudoParticle[3].GetMass() );
  pseudoParticle[5].SetTotalEnergy( pseudoParticle[3].GetTotalEnergy() );
    
  pseudoParticle[6].SetMomentum( pseudoParticle[4].GetMomentum() * (-1.0) );
  pseudoParticle[6].SetMass( pseudoParticle[4].GetMass() );
  pseudoParticle[6].SetTotalEnergy( pseudoParticle[4].GetTotalEnergy() );
  
  G4double wgt;
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  if( forwardCount > 1 )     // tempV will contain the forward particles
  {
    G4FastVector<G4ReactionProduct,256> tempV;
    tempV.Initialize( forwardCount );
    G4bool constantCrossSection = true;
    G4int tempLen = 0;
    if( currentParticle.GetSide() == 1 )
      tempV.SetElement( tempLen++, &currentParticle );
    if( targetParticle.GetSide() == 1 )
      tempV.SetElement( tempLen++, &targetParticle );
    for( i=0; i<vecLen; ++i )
    {
      if( vec[i]->GetSide() == 1 )
      {
        if( tempLen < 18 )
          tempV.SetElement( tempLen++, vec[i] );
        else
        {
          vec[i]->SetSide( -1 );
          continue;
        }
      }
    }
    if( tempLen >= 2 )
    {
      wgt = GenerateNBodyEvent( pseudoParticle[3].GetMass()/MeV,
                                constantCrossSection, tempV, tempLen );
      if( currentParticle.GetSide() == 1 )
        currentParticle.Lorentz( currentParticle, pseudoParticle[5] );
      if( targetParticle.GetSide() == 1 )
        targetParticle.Lorentz( targetParticle, pseudoParticle[5] );
      for( i=0; i<vecLen; ++i )
      {
        if( vec[i]->GetSide() == 1 )vec[i]->Lorentz( *vec[i], pseudoParticle[5] );
      }
    }
  }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  if( backwardCount > 1 )   //  tempV will contain the backward particles,
  {                         //  but not those created from the intranuclear cascade
    G4FastVector<G4ReactionProduct,256> tempV;
    tempV.Initialize( backwardCount );
    G4bool constantCrossSection = true;
    G4int tempLen = 0;
    if( currentParticle.GetSide() == -1 )
      tempV.SetElement( tempLen++, &currentParticle );
    if( targetParticle.GetSide() == -1 )
      tempV.SetElement( tempLen++, &targetParticle );
    for( i=0; i<vecLen; ++i )
    {
      if( vec[i]->GetSide() == -1 )
      {
        if( tempLen < 18 )
          tempV.SetElement( tempLen++, vec[i] );
        else
        {
          vec[i]->SetSide( -2 );
          vec[i]->SetKineticEnergy( 0.0 );
          vec[i]->SetMomentum( 0.0, 0.0, 0.0 );
          continue;
        }
      }
    }
    if( tempLen >= 2 )
    {
      wgt = GenerateNBodyEvent( pseudoParticle[4].GetMass()/MeV,
                                constantCrossSection, tempV, tempLen );
      if( currentParticle.GetSide() == -1 )
        currentParticle.Lorentz( currentParticle, pseudoParticle[6] );
      if( targetParticle.GetSide() == -1 )
        targetParticle.Lorentz( targetParticle, pseudoParticle[6] );
      for( i=0; i<vecLen; ++i )
      {
        if( vec[i]->GetSide() == -1 )vec[i]->Lorentz( *vec[i], pseudoParticle[6] );
      }
    }
  }

      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  //
  // Lorentz transformation in lab system
  //
  currentParticle.Lorentz( currentParticle, pseudoParticle[2] );
  targetParticle.Lorentz( targetParticle, pseudoParticle[2] );
  for( i=0; i<vecLen; ++i ) vec[i]->Lorentz( *vec[i], pseudoParticle[2] );

      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  //
  // sometimes the leading strange particle is lost, set it back
  //
  G4bool dum = true;
  if( leadFlag )
  {
    // leadFlag will be true
    //  iff original particle is strange AND if incident particle is strange
    //  leadFlag is set to the incident particle
    //  or
    //  target particle is strange leadFlag is set to the target particle

    if( currentParticle.GetDefinition() == leadingStrangeParticle.GetDefinition() )
      dum = false;
    else if( targetParticle.GetDefinition() == leadingStrangeParticle.GetDefinition() )
      dum = false;
    else
    {
      for( i=0; i<vecLen; ++i )
      {
        if( vec[i]->GetDefinition() == leadingStrangeParticle.GetDefinition() )
        {
          dum = false;
          break;
        }
      }
    }
    if( dum )
    {
      G4double leadMass = leadingStrangeParticle.GetMass()/MeV;
      G4double ekin;
      if( ((leadMass <  protonMass) && (targetParticle.GetMass()/MeV <  protonMass)) ||
          ((leadMass >= protonMass) && (targetParticle.GetMass()/MeV >= protonMass)) )
      {
        ekin = targetParticle.GetKineticEnergy()/GeV;
        pp1 = targetParticle.GetMomentum().mag()/MeV; // old momentum
        targetParticle.SetDefinition( leadingStrangeParticle.GetDefinition() );
        targetParticle.SetKineticEnergy( ekin*GeV );
        pp = targetParticle.GetTotalMomentum()/MeV;   // new momentum
        if( pp1 < 1.0e-3 ) {
          G4ThreeVector iso = Isotropic(pp);
          targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
        } else {
          targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
	}
        targetHasChanged = true;
      }
      else
      {
        ekin = currentParticle.GetKineticEnergy()/GeV;
        pp1 = currentParticle.GetMomentum().mag()/MeV;
        currentParticle.SetDefinition( leadingStrangeParticle.GetDefinition() );
        currentParticle.SetKineticEnergy( ekin*GeV );
        pp = currentParticle.GetTotalMomentum()/MeV;
        if( pp1 < 1.0e-3 ) {
          G4ThreeVector iso = Isotropic(pp);
          currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
        } else {
          currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
	}
        incidentHasChanged = true;
      }
    }
  }    // end of if( leadFlag )

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
  //
  //  for various reasons, the energy balance is not sufficient,
  //  check that,  energy balance, angle of final system, etc.
  //
  pseudoParticle[4].SetMass( mOriginal*GeV );
  pseudoParticle[4].SetTotalEnergy( etOriginal*GeV );
  pseudoParticle[4].SetMomentum( 0.0, 0.0, pOriginal*GeV );
    
  const G4ParticleDefinition* aOrgDef = modifiedOriginal.GetDefinition();
  G4int diff = 0;
  if(aOrgDef == G4Proton::Proton() || aOrgDef == G4Neutron::Neutron() )  diff = 1;
  if(numberofFinalStateNucleons == 1) diff = 0;
  pseudoParticle[5].SetMomentum( 0.0, 0.0, 0.0 );
  pseudoParticle[5].SetMass( protonMass*(numberofFinalStateNucleons-diff)*MeV);
  pseudoParticle[5].SetTotalEnergy( protonMass*(numberofFinalStateNucleons-diff)*MeV);
    
  G4double theoreticalKinetic =
    pseudoParticle[4].GetTotalEnergy()/GeV + pseudoParticle[5].GetTotalEnergy()/GeV;
    
  pseudoParticle[6] = pseudoParticle[4] + pseudoParticle[5];
  pseudoParticle[4].Lorentz( pseudoParticle[4], pseudoParticle[6] );
  pseudoParticle[5].Lorentz( pseudoParticle[5], pseudoParticle[6] );

  if( vecLen < 16 )
  {
    G4ReactionProduct tempR[130];
    tempR[0] = currentParticle;
    tempR[1] = targetParticle;
    for( i=0; i<vecLen; ++i )tempR[i+2] = *vec[i];

    G4FastVector<G4ReactionProduct,256> tempV;
    tempV.Initialize( vecLen+2 );
    G4bool constantCrossSection = true;
    G4int tempLen = 0;
    for( i=0; i<vecLen+2; ++i )tempV.SetElement( tempLen++, &tempR[i] );

    if( tempLen >= 2 )
    {
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      wgt = GenerateNBodyEvent( pseudoParticle[4].GetTotalEnergy()/MeV +
                                pseudoParticle[5].GetTotalEnergy()/MeV,
                                constantCrossSection, tempV, tempLen );
      if (wgt == -1) {
        G4double Qvalue = 0;
        for (i = 0; i < tempLen; i++) Qvalue += tempV[i]->GetMass();
        wgt = GenerateNBodyEvent( Qvalue/MeV,
                                  constantCrossSection, tempV, tempLen );
      }
      theoreticalKinetic = 0.0;
      for( i=0; i<vecLen+2; ++i )
      {
        pseudoParticle[7].SetMomentum( tempV[i]->GetMomentum() );
        pseudoParticle[7].SetMass( tempV[i]->GetMass() );
        pseudoParticle[7].SetTotalEnergy( tempV[i]->GetTotalEnergy() );
        pseudoParticle[7].Lorentz( pseudoParticle[7], pseudoParticle[5] );
        theoreticalKinetic += pseudoParticle[7].GetKineticEnergy()/GeV;
      }
    }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  }
  else
  {
    theoreticalKinetic -=
      ( currentParticle.GetMass()/GeV + targetParticle.GetMass()/GeV );
    for( i=0; i<vecLen; ++i )theoreticalKinetic -= vec[i]->GetMass()/GeV;
  }
  G4double simulatedKinetic =
    currentParticle.GetKineticEnergy()/GeV + targetParticle.GetKineticEnergy()/GeV;
  for( i=0; i<vecLen; ++i )simulatedKinetic += vec[i]->GetKineticEnergy()/GeV;

  // make sure that kinetic energies are correct
  // the backward nucleon cluster is not produced within proper kinematics!!!
    
  if( simulatedKinetic != 0.0 )
  {
    wgt = (theoreticalKinetic)/simulatedKinetic;
    currentParticle.SetKineticEnergy( wgt*currentParticle.GetKineticEnergy() );
    pp = currentParticle.GetTotalMomentum()/MeV;
    pp1 = currentParticle.GetMomentum().mag()/MeV;
    if( pp1 < 0.001*MeV ) {
      G4ThreeVector iso = Isotropic(pp);
      currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
    } else {
      currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
    }

    targetParticle.SetKineticEnergy( wgt*targetParticle.GetKineticEnergy() );
    pp = targetParticle.GetTotalMomentum()/MeV;
    pp1 = targetParticle.GetMomentum().mag()/MeV;
    if( pp1 < 0.001*MeV ) {
      G4ThreeVector iso = Isotropic(pp);
      targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
    } else {
      targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
    }

    for( i=0; i<vecLen; ++i )
    {
      vec[i]->SetKineticEnergy( wgt*vec[i]->GetKineticEnergy() );
      pp = vec[i]->GetTotalMomentum()/MeV;
      pp1 = vec[i]->GetMomentum().mag()/MeV;
      if( pp1 < 0.001 ) {
        G4ThreeVector iso = Isotropic(pp);
        vec[i]->SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        vec[i]->SetMomentum( vec[i]->GetMomentum() * (pp/pp1) );
      }
    }
  }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);

  Rotate( numberofFinalStateNucleons, pseudoParticle[4].GetMomentum(),
          modifiedOriginal, originalIncident, targetNucleus,
          currentParticle, targetParticle, vec, vecLen );

  //  Add black track particles
  //  the total number of particles produced is restricted to 198
  //  this may have influence on very high energies

  if( atomicWeight >= 1.5 )
  {
    // npnb is number of proton/neutron black track particles
    // ndta is the number of deuterons, tritons, and alphas produced
    // epnb is the kinetic energy available for proton/neutron black track 
    //   particles
    // edta is the kinetic energy available for deuteron/triton/alpha 
    //   particles

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

    const G4double pnCutOff = 0.001;     // GeV
    const G4double dtaCutOff = 0.001;    // GeV
    //    const G4double kineticMinimum = 1.e-6;
    //    const G4double kineticFactor = -0.005;
      
    //    G4double sprob = 0.0;   // sprob = probability of self-absorption in 
                            // heavy molecules
    // Not currently used (DHW 9 June 2008)  const G4double ekIncident = originalIncident->GetKineticEnergy()/GeV;
    //    if( ekIncident >= 5.0 )sprob = std::min( 1.0, 0.6*std::log(ekIncident-4.0) );
      
    if( epnb >= pnCutOff )
    {
      npnb = G4Poisson((1.5+1.25*numberofFinalStateNucleons)*epnb/(epnb+edta));
      if( numberofFinalStateNucleons + npnb > atomicWeight )
        npnb = G4int(atomicWeight - numberofFinalStateNucleons);
      npnb = std::min( npnb, 127-vecLen );
    }
    if( edta >= dtaCutOff )
    {
      ndta = G4Poisson( (1.5+1.25*numberofFinalStateNucleons)*edta/(epnb+edta) );
      ndta = std::min( ndta, 127-vecLen );
    }
    if (npnb == 0 && ndta == 0) npnb = 1;

    // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);

    AddBlackTrackParticles(epnb, npnb, edta, ndta, modifiedOriginal, 
                           PinNucleus, NinNucleus, targetNucleus,
                           vec, vecLen );
    // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
  }

  //if( centerofmassEnergy <= (4.0+G4UniformRand()) )
  //  MomentumCheck( modifiedOriginal, currentParticle, targetParticle, vec, vecLen );
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
