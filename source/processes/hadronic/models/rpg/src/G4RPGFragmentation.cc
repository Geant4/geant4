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
// $Id: G4RPGFragmentation.cc,v 1.3 2007/12/06 01:13:14 dennis Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
 
#include "G4RPGFragmentation.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include <iostream>
#include "G4HadReentrentException.hh"
#include <signal.h>


G4RPGFragmentation::G4RPGFragmentation()
  : G4RPGReaction() {}


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
  // Derived from H. Fesefeldt's original FORTRAN code GENXPT
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
    if(vecLen == 0) return false;

    G4ParticleDefinition *aPiMinus = G4PionMinus::PionMinus();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aPiPlus = G4PionPlus::PionPlus();
    G4ParticleDefinition *aPiZero = G4PionZero::PionZero();

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
    //
    //  randomize the order of the secondary particles
    //  note that the current and target particles are not affected
    //
    for( i=0; i<vecLen; ++i )
    {
      G4int itemp = G4int( G4UniformRand()*vecLen );
      G4ReactionProduct pTemp = *vec[itemp];
      *vec[itemp] = *vec[i];
      *vec[i] = pTemp;
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    }

    if( currentMass == 0.0 && targetMass == 0.0 )
    {
      // Target and projectile have annihilated.  Replace them with the first 
      // two secondaries in the list.  Current particle KE is maintained.
 
      G4double ek = currentParticle.GetKineticEnergy();
      G4ThreeVector m = currentParticle.GetMomentum();
      currentParticle = *vec[0];
      targetParticle = *vec[1];
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
      currentParticle.SetMomentum( m );
      veryForward = true;
    }
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    const G4double protonMass = aProton->GetPDGMass()/MeV;

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
          0.312+0.200*std::log(std::log(centerofmassEnergy*centerofmassEnergy))+
          std::pow(centerofmassEnergy*centerofmassEnergy,1.5)/6000.0 );
    
    G4double freeEnergy = centerofmassEnergy-currentMass-targetMass;
    G4double forwardEnergy = freeEnergy/2.;
    G4int forwardCount = 1;         // number of particles in forward hemisphere
    
    G4double backwardEnergy = freeEnergy/2.;
    G4int backwardCount = 1;        // number of particles in backward hemisphere

    if(veryForward)
    {
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

    for( i=0; i<vecLen; ++i )
    {
      if( vec[i]->GetSide() == -1 )
      {
        ++backwardCount;
        backwardEnergy -= vec[i]->GetMass()/GeV;
      } else {
        ++forwardCount;
        forwardEnergy -= vec[i]->GetMass()/GeV;
      }
    }
    //
    //  Add particles from intranuclear cascade.
    //  nuclearExcitationCount = number of new secondaries produced by nuclear excitation
    //  extraCount = number of nucleons within these new secondaries
    //
    G4double xtarg;
    if( centerofmassEnergy < (2.0+G4UniformRand()) )
      xtarg = afc * (std::pow(atomicWeight,0.33)-1.0) * (2.0*backwardCount+vecLen+2)/2.0;
    else
      xtarg = afc * (std::pow(atomicWeight,0.33)-1.0) * (2.0*backwardCount);
    if( xtarg <= 0.0 )xtarg = 0.01;
    G4int nuclearExcitationCount = G4Poisson( xtarg );
    if(atomicWeight<1.0001) nuclearExcitationCount = 0;
    G4int extraNucleonCount = 0;
    G4double extraNucleonMass = 0.0;
    if( nuclearExcitationCount > 0 )
    {
      const G4double nucsup[] = { 1.00, 0.7, 0.5, 0.4, 0.35, 0.3 };
      const G4double psup[] = { 3., 6., 20., 50., 100., 1000. };
      G4int momentumBin = 0;
      while( (momentumBin < 6) &&
             (modifiedOriginal.GetTotalMomentum()/GeV > psup[momentumBin]) )
        ++momentumBin;
      momentumBin = std::min( 5, momentumBin );
      //
      //  NOTE: in GENXPT, these new particles were given negative codes
      //        here I use  NewlyAdded = true  instead
      //

      for( i=0; i<nuclearExcitationCount; ++i )
      {
        G4ReactionProduct * pVec = new G4ReactionProduct();
	if( G4UniformRand() < nucsup[momentumBin] )
        {
          if( G4UniformRand() > 1.0-atomicNumber/atomicWeight )
            pVec->SetDefinition( aProton );
          else
            pVec->SetDefinition( aNeutron );
          pVec->SetSide( -2 );                // -2 means backside nucleon
          ++extraNucleonCount;
          backwardEnergy += pVec->GetMass()/GeV;
          extraNucleonMass += pVec->GetMass()/GeV;
        }
        else
        {
          G4double ran = G4UniformRand();
          if( ran < 0.3181 )
            pVec->SetDefinition( aPiPlus );
          else if( ran < 0.6819 )
            pVec->SetDefinition( aPiZero );
          else
            pVec->SetDefinition( aPiMinus );
          pVec->SetSide( -1 );                // backside particle, but not a nucleon
        }
        pVec->SetNewlyAdded( true );                // true is the same as IPA(i)<0
        vec.SetElement( vecLen++, pVec );    
        // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
        backwardEnergy -= pVec->GetMass()/GeV;
        ++backwardCount;
      }
    }
    //
    //  assume conservation of kinetic energy in forward & backward hemispheres
    //
    G4int is, iskip;
    while( forwardEnergy <= 0.0 )  // must eliminate a particle from the forward side
    {
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      iskip = G4int(G4UniformRand()*forwardCount) + 1; // 1 <= iskip <= forwardCount
      is = 0;
      G4int forwardParticlesLeft = 0;
      for( i=(vecLen-1); i>=0; --i )
      {
        if( vec[i]->GetSide() == 1 && vec[i]->GetMayBeKilled())
        {
          forwardParticlesLeft = 1;  
          if( ++is == iskip )
          { 
            forwardEnergy += vec[i]->GetMass()/GeV;
            for( G4int j=i; j<(vecLen-1); j++ )*vec[j] = *vec[j+1];    // shift up
            --forwardCount;
            delete vec[vecLen-1];
            if( --vecLen == 0 )return false;  // all the secondaries have been eliminated
            break;  // --+
          }         //   |
        }           //   |
      }             // break goes down to here
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      if( forwardParticlesLeft == 0 )
      {
        G4int iremove = -1;
        for (G4int i = 0; i < vecLen; i++) {
          if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
            iremove = i;
            break;
          }
        }
        if (iremove == -1) {
          for (G4int i = 0; i < vecLen; i++) {
            if (vec[i]->GetDefinition()->GetParticleSubType() == "kaon") {
              iremove = i;
              break;
            }
          }
        }
        if (iremove == -1) iremove = 0;

        forwardEnergy += vec[iremove]->GetMass()/GeV;
        if (vec[iremove]->GetSide() > 0) --forwardCount;
 
        for (G4int i = iremove; i < vecLen-1; i++) *vec[i] = *vec[i+1];
        delete vec[vecLen-1];
        vecLen--;
        if (vecLen == 0) return false;  // all secondaries have been eliminated
        break;
      }
    } // while

      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    while( backwardEnergy <= 0.0 )  // must eliminate a particle from the backward side
    {
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      iskip = G4int(G4UniformRand()*backwardCount) + 1; // 1 <= iskip <= backwardCount
      is = 0;
      G4int backwardParticlesLeft = 0;
      for( i=(vecLen-1); i>=0; --i )
      {
        if( vec[i]->GetSide() < 0 && vec[i]->GetMayBeKilled())
	{
          backwardParticlesLeft = 1;
          if( ++is == iskip )        // eliminate the i'th particle
          {
            if( vec[i]->GetSide() == -2 )
            {
              --extraNucleonCount;
              extraNucleonMass -= vec[i]->GetMass()/GeV;
              backwardEnergy -= vec[i]->GetTotalEnergy()/GeV;
            }
            backwardEnergy += vec[i]->GetTotalEnergy()/GeV;
            for( G4int j=i; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];   // shift up
            --backwardCount;
            delete vec[vecLen-1];
            if( --vecLen == 0 )return false;  // all the secondaries have been eliminated
            break;
          }
        }
      }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      if( backwardParticlesLeft == 0 )
      {
        G4int iremove = -1;
        for (G4int i = 0; i < vecLen; i++) {
          if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
            iremove = i;
            break;
          }
        }
        if (iremove == -1) {
          for (G4int i = 0; i < vecLen; i++) {
            if (vec[i]->GetDefinition()->GetParticleSubType() == "kaon") {
              iremove = i;
              break;
            }
          }
        }
        if (iremove == -1) iremove = 0;
 
        backwardEnergy += vec[iremove]->GetMass()/GeV;
        if (vec[iremove]->GetSide() > 0) --backwardCount;
 
        for (G4int i = iremove; i < vecLen-1; i++) *vec[i] = *vec[i+1];
        delete vec[vecLen-1];
        vecLen--;
        if (vecLen == 0) return false;  // all secondaries have been eliminated
        break;
      }
    } // while

      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    //
    //  define initial state vectors for Lorentz transformations
    //  the pseudoParticles have non-standard masses, hence the "pseudo"
    //
    G4ReactionProduct pseudoParticle[8];
    for( i=0; i<8; ++i )pseudoParticle[i].SetZero();
    
    pseudoParticle[0].SetMass( mOriginal*GeV );
    pseudoParticle[0].SetMomentum( 0.0, 0.0, pOriginal*GeV );
    pseudoParticle[0].SetTotalEnergy(
     std::sqrt( pOriginal*pOriginal + mOriginal*mOriginal )*GeV );
    
    pseudoParticle[1].SetMass( protonMass*MeV );        // this could be targetMass
    pseudoParticle[1].SetTotalEnergy( protonMass*MeV );
    
    pseudoParticle[3].SetMass( protonMass*(1+extraNucleonCount)*MeV );
    pseudoParticle[3].SetTotalEnergy( protonMass*(1+extraNucleonCount)*MeV );
    
    pseudoParticle[2] = pseudoParticle[0] + pseudoParticle[1];
    pseudoParticle[3] = pseudoParticle[3] + pseudoParticle[0];
    
    pseudoParticle[0].Lorentz( pseudoParticle[0], pseudoParticle[2] );
    pseudoParticle[1].Lorentz( pseudoParticle[1], pseudoParticle[2] );
    
    //
    //  main loop for 4-momentum generation
    //  see Pitha-report (Aachen) for a detailed description of the method
    //
    G4double aspar, pt, et, x, pp, pp1, wgt;
    G4int    innerCounter, outerCounter;
    G4bool   eliminateThisParticle, resetEnergies, constantCrossSection;
    
    G4double forwardKinetic = 0.0, backwardKinetic = 0.0;
    //
    // process the secondary particles in reverse order
    // the incident particle is Done after the secondaries
    // nucleons, including the target, in the backward hemisphere are also Done later
    //
    G4int backwardNucleonCount = 0;       // number of nucleons in backward hemisphere
    G4double totalEnergy, kineticEnergy, vecMass;

    for( i=(vecLen-1); i>=0; --i )
    {
      G4double phi = G4UniformRand()*twopi;
      if( vec[i]->GetNewlyAdded() )           // added from intranuclear cascade
      {
        if( vec[i]->GetSide() == -2 )         //  is a nucleon
        {
          if( backwardNucleonCount < 18 )
          {
            if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
              for(G4int i=0; i<vecLen; i++) delete vec[i];
   	      vecLen = 0;
              throw G4HadReentrentException(__FILE__, __LINE__,
	      "G4RPGFragmentation::ReactionStage : a pion has been counted as a backward nucleon");
            }
            vec[i]->SetSide( -3 );
            ++backwardNucleonCount;
            continue;
          }
        }
      }
      //
      //  set pt and phi values, they are changed somewhat in the iteration loop
      //  set mass parameter for lambda fragmentation model
      //
      vecMass = vec[i]->GetMass()/GeV;
      G4double ran = -std::log(1.0-G4UniformRand())/3.5;
      if( vec[i]->GetSide() == -2 )   // backward nucleon
      {
        if (vec[i]->GetDefinition()->GetParticleSubType() == "kaon" ||
            vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
          aspar = 0.75;
          pt = std::sqrt( std::pow( ran, 1.7 ) );
        } else {                          // vec[i] must be a proton, neutron,
          aspar = 0.20;                   //  lambda, sigma, xsi, or ion
          pt = std::sqrt( std::pow( ran, 1.2 ) );
        }

      } else {                          // not a backward nucleon

        if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
          aspar = 0.75;
          pt = std::sqrt( std::pow( ran, 1.7 ) );
        } else if (vec[i]->GetDefinition()->GetParticleSubType() == "kaon") {
          aspar = 0.70;
          pt = std::sqrt( std::pow( ran, 1.7 ) );
        } else {                        // vec[i] must be a proton, neutron, 
          aspar = 0.65;                 //  lambda, sigma, xsi, or ion
          pt = std::sqrt( std::pow( ran, 1.5 ) );
        }
      }
      pt = std::max( 0.001, pt );
      vec[i]->SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );
      if( vec[i]->GetSide() > 0 )
        et = pseudoParticle[0].GetTotalEnergy()/GeV;
      else
        et = pseudoParticle[1].GetTotalEnergy()/GeV;

      //
      //   start of outer iteration loop
      //
      outerCounter = 0;
      eliminateThisParticle = true;
      resetEnergies = true;
      dndl[0] = 0.0;

      while( ++outerCounter < 3 )
      {
        FragmentationIntegral(pt, et, aspar, vecMass);

        innerCounter = 0;
        vec[i]->SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );
        //
        //   start of inner iteration loop
        //
        while( ++innerCounter < 7 )
        {
          ran = G4UniformRand()*dndl[19];
          l = 1;
          while( ( ran > dndl[l] ) && ( l < 19 ) ) l++;
          x = (G4double(l-1) + G4UniformRand())/19.;
          if( vec[i]->GetSide() < 0 )x *= -1.;
          vec[i]->SetMomentum( x*et*GeV );              // set the z-momentum
          totalEnergy = std::sqrt( x*et*x*et + pt*pt + vecMass*vecMass );
          vec[i]->SetTotalEnergy( totalEnergy*GeV );
          kineticEnergy = vec[i]->GetKineticEnergy()/GeV;
          if( vec[i]->GetSide() > 0 )                            // forward side
          {
            if( (forwardKinetic+kineticEnergy) < 0.95*forwardEnergy )
            {
              pseudoParticle[4] = pseudoParticle[4] + (*vec[i]);
              forwardKinetic += kineticEnergy;
              outerCounter = 2;                     // leave outer loop
              eliminateThisParticle = false;        // don't eliminate this particle
              resetEnergies = false;
              break;                                // leave inner loop
            }
            if( innerCounter > 5 )break;           // leave inner loop
            if( backwardEnergy >= vecMass )  // switch sides
            {
              vec[i]->SetSide( -1 );
              forwardEnergy += vecMass;
              backwardEnergy -= vecMass;
              ++backwardCount;
            }
          } else {                                                 // backward side
           if( extraNucleonCount > 19 )   // commented out to duplicate ?bug? in GENXPT
             x = 0.999;
           G4double xxx = 0.95+0.05*extraNucleonCount/20.0;
           if( (backwardKinetic+kineticEnergy) < xxx*backwardEnergy )
            {
              pseudoParticle[5] = pseudoParticle[5] + (*vec[i]);
              backwardKinetic += kineticEnergy;
              outerCounter = 2;                    // leave outer loop
              eliminateThisParticle = false;       // don't eliminate this particle
              resetEnergies = false;
              break;                               // leave inner loop
            }
            if( innerCounter > 5 )break;           // leave inner loop
            if( forwardEnergy >= vecMass )  // switch sides
            {
              vec[i]->SetSide( 1 );
              forwardEnergy -= vecMass;
              backwardEnergy += vecMass;
              backwardCount--;
            }
          }
          G4ThreeVector momentum = vec[i]->GetMomentum();
          vec[i]->SetMomentum( momentum.x() * 0.9, momentum.y() * 0.9 );
          pt *= 0.9;
          dndl[19] *= 0.9;
        }                       // closes inner loop
        if( resetEnergies ) {
          //  If we get to here, the inner loop has been done 6 times.
          //  Reset the kinetic energies of previously done particles, if 
	  //  they are lighter than protons and in the forward hemisphere,
          //  then continue with outer loop.
          //
          forwardKinetic = 0.0;
          backwardKinetic = 0.0;
          pseudoParticle[4].SetZero();
          pseudoParticle[5].SetZero();
          for( l=i+1; l<vecLen; ++l ) {
            if (vec[l]->GetSide() > 0 ||
                vec[l]->GetDefinition()->GetParticleSubType() == "kaon" ||
                vec[l]->GetDefinition()->GetParticleSubType() == "pi") {

              G4double tempMass = vec[l]->GetMass()/MeV;
              totalEnergy = 0.95*vec[l]->GetTotalEnergy()/MeV + 0.05*tempMass;
              totalEnergy = std::max( tempMass, totalEnergy );
              vec[l]->SetTotalEnergy( totalEnergy*MeV );
              pp = std::sqrt( std::abs( totalEnergy*totalEnergy - tempMass*tempMass ) );
              pp1 = vec[l]->GetMomentum().mag()/MeV;
              if( pp1 < 1.0e-6*GeV ) {
                G4ThreeVector iso = Isotropic(pp);
                vec[l]->SetMomentum( iso.x(), iso.y(), iso.z() );
              } else {
                vec[l]->SetMomentum( vec[l]->GetMomentum() * (pp/pp1) );
              }
              G4double px = vec[l]->GetMomentum().x()/MeV;
              G4double py = vec[l]->GetMomentum().y()/MeV;
              pt = std::max( 1.0, std::sqrt( px*px + py*py ) )/GeV;
              if( vec[l]->GetSide() > 0 )
              {
                forwardKinetic += vec[l]->GetKineticEnergy()/GeV;
                pseudoParticle[4] = pseudoParticle[4] + (*vec[l]);
              } else {
                backwardKinetic += vec[l]->GetKineticEnergy()/GeV;
                pseudoParticle[5] = pseudoParticle[5] + (*vec[l]);
              }
            } // if pi, K or forward
          } // for l
        } // if resetEnergies
      } // closes outer loop
              
      if( eliminateThisParticle && vec[i]->GetMayBeKilled())  // not enough energy, eliminate this particle
      {
        if( vec[i]->GetSide() > 0 )
        {
          --forwardCount;
          forwardEnergy += vecMass;
        } else {
          if( vec[i]->GetSide() == -2 )
          {
            --extraNucleonCount;
            extraNucleonMass -= vecMass;
            backwardEnergy -= vecMass;
          }
          --backwardCount;
          backwardEnergy += vecMass;
        }
        for( G4int j=i; j<(vecLen-1); ++j )*vec[j] = *vec[j+1];    // shift up
        G4ReactionProduct *temp = vec[vecLen-1];
        delete temp;
        // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
        if( --vecLen == 0 )return false;  // all the secondaries have been eliminated
      }
    }   // closes main for loop

    //
    //  for the incident particle:  it was placed in the forward hemisphere
    //   set pt and phi values, they are changed somewhat in the iteration loop
    //   set mass parameter for lambda fragmentation model
    //
    G4double phi = G4UniformRand()*twopi;
    G4double ran = -std::log(1.0-G4UniformRand());
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

    currentParticle.SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );
    et = pseudoParticle[0].GetTotalEnergy()/GeV;
    dndl[0] = 0.0;
    vecMass = currentParticle.GetMass()/GeV;

    FragmentationIntegral(pt, et, aspar, vecMass);

    ran = G4UniformRand()*dndl[19];
    l = 1;
    while( ( ran > dndl[l] ) && ( l < 19 ) ) l++;
    x = (G4double(l-1) + G4UniformRand())/19.;
    currentParticle.SetMomentum( x*et*GeV );                 // set the z-momentum
    if( forwardEnergy < forwardKinetic )
      totalEnergy = vecMass + 0.04*std::fabs(normal());
    else
      totalEnergy = vecMass + forwardEnergy - forwardKinetic;
    currentParticle.SetTotalEnergy( totalEnergy*GeV );
    pp = std::sqrt( std::abs( totalEnergy*totalEnergy - vecMass*vecMass ) )*GeV;
    pp1 = currentParticle.GetMomentum().mag()/MeV;

    if( pp1 < 1.0e-6*GeV ) {
      G4ThreeVector iso = Isotropic(pp);
      currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
    } else {
      currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
    }
    pseudoParticle[4] = pseudoParticle[4] + currentParticle;

    //
    // Current particle now finished
    //
    // Begin target particle
    //

    if( backwardNucleonCount < 18 )
    {
      targetParticle.SetSide( -3 );
      ++backwardNucleonCount;
    }
    else
    {
      //  set pt and phi values, they are changed somewhat in the iteration loop
      //  set mass parameter for lambda fragmentation model
      //
      vecMass = targetParticle.GetMass()/GeV;
      ran = -std::log(1.0-G4UniformRand());
      aspar = 0.40;
      pt = std::max( 0.001, std::sqrt( std::pow( ran/4.0, 1.2 ) ) );
      targetParticle.SetMomentum( pt*std::cos(phi)*GeV, pt*std::sin(phi)*GeV );
      et = pseudoParticle[1].GetTotalEnergy()/GeV;
      outerCounter = 0;
      eliminateThisParticle = true;     // should never eliminate the target particle
      resetEnergies = true;
      dndl[0] = 0.0;

      while( ++outerCounter < 3 )     // start of outer iteration loop
      {
        FragmentationIntegral(pt, et, aspar, vecMass);

        innerCounter = 0;
        while( ++innerCounter < 7 )    // start of inner iteration loop
        {
          ran = G4UniformRand()*dndl[19];
          l = 1;
          while( ( ran > dndl[l] ) && ( l < 19 ) ) l++;
          x = (G4double(l-1) + G4UniformRand())/19.;
          if( targetParticle.GetSide() < 0 )x *= -1.;
          targetParticle.SetMomentum( x*et*GeV );                // set the z-momentum
          totalEnergy = std::sqrt( x*et*x*et + pt*pt + vecMass*vecMass );
          targetParticle.SetTotalEnergy( totalEnergy*GeV );
          if( targetParticle.GetSide() < 0 )
          {
            if( extraNucleonCount > 19 )x=0.999;
            G4double xxx = 0.95+0.05*extraNucleonCount/20.0;
            if( (backwardKinetic+totalEnergy-vecMass) < xxx*backwardEnergy )
            {
              pseudoParticle[5] = pseudoParticle[5] + targetParticle;
              backwardKinetic += totalEnergy - vecMass;
	      //              pseudoParticle[6] = pseudoParticle[4] + pseudoParticle[5];
	      //              pseudoParticle[6].SetMomentum( 0.0 );                      // set z-momentum
              outerCounter = 2;                    // leave outer loop
              eliminateThisParticle = false;       // don't eliminate this particle
              resetEnergies = false;
              break;                               // leave inner loop
            }
            if( innerCounter > 5 )break;           // leave inner loop
            if( forwardEnergy >= vecMass )  // switch sides
            {
              targetParticle.SetSide( 1 );
              forwardEnergy -= vecMass;
              backwardEnergy += vecMass;
              --backwardCount;
            }
            G4ThreeVector momentum = targetParticle.GetMomentum();
            targetParticle.SetMomentum( momentum.x() * 0.9, momentum.y() * 0.9 );
            pt *= 0.9;
            dndl[19] *= 0.9;
          }
          else                    // target has gone to forward side
          {
            if( forwardEnergy < forwardKinetic )
              totalEnergy = vecMass + 0.04*std::fabs(normal());
            else
              totalEnergy = vecMass + forwardEnergy - forwardKinetic;
            targetParticle.SetTotalEnergy( totalEnergy*GeV );
            pp = std::sqrt( std::abs( totalEnergy*totalEnergy - vecMass*vecMass ) )*GeV;
            pp1 = targetParticle.GetMomentum().mag()/MeV;
            if( pp1 < 1.0e-6*GeV ) {
              G4ThreeVector iso = Isotropic(pp);
              targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
            } else {
              targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
	    }

            pseudoParticle[4] = pseudoParticle[4] + targetParticle;
            outerCounter = 2;                    // leave outer loop
            eliminateThisParticle = false;       // don't eliminate this particle
            resetEnergies = false;
            break;                               // leave inner loop
          }
        }    // closes inner loop

        if( resetEnergies ) {
          //  If we get to here, the inner loop has been done 6 times.
          //  Reset the kinetic energies of previously done particles, 
	  //  if they are lighter than protons and in the forward hemisphere,
          //  then continue with outer loop.
          
          forwardKinetic = backwardKinetic = 0.0;
          pseudoParticle[4].SetZero();
          pseudoParticle[5].SetZero();
          for( l=0; l<vecLen; ++l ) {
            if (vec[l]->GetSide() > 0 ||
                vec[l]->GetDefinition()->GetParticleSubType() == "kaon" ||
                vec[l]->GetDefinition()->GetParticleSubType() == "pi") {
              G4double tempMass = vec[l]->GetMass()/GeV;
              totalEnergy =
                std::max( tempMass, 0.95*vec[l]->GetTotalEnergy()/GeV + 0.05*tempMass );
              vec[l]->SetTotalEnergy( totalEnergy*GeV );
              pp = std::sqrt( std::abs( totalEnergy*totalEnergy - tempMass*tempMass ) )*GeV;
              pp1 = vec[l]->GetMomentum().mag()/MeV;
              if( pp1 < 1.0e-6*GeV ) {
                G4ThreeVector iso = Isotropic(pp);
                vec[l]->SetMomentum( iso.x(), iso.y(), iso.z() );
              } else {
                vec[l]->SetMomentum( vec[l]->GetMomentum() * (pp/pp1) );
	      }
              pt = std::max( 0.001*GeV, std::sqrt( sqr(vec[l]->GetMomentum().x()/MeV) +
                                         sqr(vec[l]->GetMomentum().y()/MeV) ) )/GeV;
              if( vec[l]->GetSide() > 0)
              {
                forwardKinetic += vec[l]->GetKineticEnergy()/GeV;
                pseudoParticle[4] = pseudoParticle[4] + (*vec[l]);
              } else {
                backwardKinetic += vec[l]->GetKineticEnergy()/GeV;
                pseudoParticle[5] = pseudoParticle[5] + (*vec[l]);
              }
            } // if pi, K or forward
          } // for l
        } // if (resetEnergies)
      } // closes outer loop

//      if( eliminateThisParticle )  // not enough energy, eliminate target
//      {
//        G4cerr << "Warning: eliminating target particle" << G4endl;
//        exit( EXIT_FAILURE );
//      }
    }
    //
    // Target particle finished.
    //
    // Now produce backward nucleons with a cluster model
    //
    pseudoParticle[6].Lorentz( pseudoParticle[3], pseudoParticle[2] );
    pseudoParticle[6] = pseudoParticle[6] - pseudoParticle[4];
    pseudoParticle[6] = pseudoParticle[6] - pseudoParticle[5];
    if( backwardNucleonCount == 1 )  // target particle is the only backward nucleon
    {
      G4double ekin =
        std::min( backwardEnergy-backwardKinetic, centerofmassEnergy/2.0-protonMass/GeV );

      if( ekin < 0.04 )ekin = 0.04 * std::fabs( normal() );
      vecMass = targetParticle.GetMass()/GeV;
      totalEnergy = ekin+vecMass;
      targetParticle.SetTotalEnergy( totalEnergy*GeV );
      pp = std::sqrt( std::abs( totalEnergy*totalEnergy - vecMass*vecMass ) )*GeV;
      pp1 = pseudoParticle[6].GetMomentum().mag()/MeV;
      if( pp1 < 1.0e-6*GeV ) { 
        G4ThreeVector iso = Isotropic(pp);
        targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        targetParticle.SetMomentum( pseudoParticle[6].GetMomentum() * (pp/pp1) );
      }
      pseudoParticle[5] = pseudoParticle[5] + targetParticle;
    }
    else if (backwardNucleonCount > 1)
    {
      const G4double cpar[] = { 1.60, 1.35, 1.15, 1.10 };
      const G4double gpar[] = { 2.60, 1.80, 1.30, 1.20 };

      G4int tempCount = 5;
      if (backwardNucleonCount < 5) tempCount = backwardNucleonCount;
      tempCount -= 2;

      G4double rmb = 0.;
      if( targetParticle.GetSide() == -3 ) rmb += targetParticle.GetMass()/GeV;
      for( i=0; i<vecLen; ++i )
      {
        if( vec[i]->GetSide() == -3 ) rmb += vec[i]->GetMass()/GeV;
      }
      rmb += std::pow(-std::log(1.0-G4UniformRand()),1./cpar[tempCount]) / gpar[tempCount];
      totalEnergy = pseudoParticle[6].GetTotalEnergy()/GeV;
      vecMass = std::min( rmb, totalEnergy );
      pseudoParticle[6].SetMass( vecMass*GeV );
      pp = std::sqrt( std::abs( totalEnergy*totalEnergy - vecMass*vecMass ) )*GeV;
      pp1 = pseudoParticle[6].GetMomentum().mag()/MeV;
      if( pp1 < 1.0e-6*GeV ) {
        G4ThreeVector iso = Isotropic(pp);
        pseudoParticle[6].SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        pseudoParticle[6].SetMomentum( pseudoParticle[6].GetMomentum() * (-pp/pp1) );
      }
      G4FastVector<G4ReactionProduct,256> tempV;  // tempV contains the backward nucleons
      tempV.Initialize( backwardNucleonCount );
      G4int tempLen = 0;
      if( targetParticle.GetSide() == -3 )tempV.SetElement( tempLen++, &targetParticle );
      for( i=0; i<vecLen; ++i )
      {
        if( vec[i]->GetSide() == -3 )tempV.SetElement( tempLen++, vec[i] );
      }
      if( tempLen != backwardNucleonCount )
      {
        G4cerr << "tempLen is not the same as backwardNucleonCount" << G4endl;
        G4cerr << "tempLen = " << tempLen;
        G4cerr << ", backwardNucleonCount = " << backwardNucleonCount << G4endl;
        G4cerr << "targetParticle side = " << targetParticle.GetSide() << G4endl;
        G4cerr << "currentParticle side = " << currentParticle.GetSide() << G4endl;
        for( i=0; i<vecLen; ++i )
          G4cerr << "particle #" << i << " side = " << vec[i]->GetSide() << G4endl;
        exit( EXIT_FAILURE );
      }
      constantCrossSection = true;
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      if( tempLen >= 2 )
      {
        wgt = GenerateNBodyEvent(
         pseudoParticle[6].GetMass(), constantCrossSection, tempV, tempLen );
         // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
        if( targetParticle.GetSide() == -3 )
        {
          targetParticle.Lorentz( targetParticle, pseudoParticle[6] );
          // tempV contains the real stuff
          pseudoParticle[5] = pseudoParticle[5] + targetParticle;
        }
        for( i=0; i<vecLen; ++i )
        {
          if( vec[i]->GetSide() == -3 )
          {
            vec[i]->Lorentz( *vec[i], pseudoParticle[6] );
            pseudoParticle[5] = pseudoParticle[5] + (*vec[i]);
          }
        }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      }
    }
    else return false;

    //
    //  Lorentz transformation in lab system
    //
    if( vecLen == 0 )return false;  // all the secondaries have been eliminated
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    
    currentParticle.Lorentz( currentParticle, pseudoParticle[1] );
    targetParticle.Lorentz( targetParticle, pseudoParticle[1] );    
    for( i=0; i<vecLen; ++i ) vec[i]->Lorentz( *vec[i], pseudoParticle[1] );

    // leadFlag will be true if original particle and incident particle are
    // both strange, in which case the incident particle becomes the leading 
    // particle. 
    // leadFlag will also be true if the target particle is strange, but the
    // incident particle is not, in which case the target particle becomes the 
    // leading particle.

    G4bool leadingStrangeParticleHasChanged = true;
    if( leadFlag )
    {
      if( currentParticle.GetDefinition() == leadingStrangeParticle.GetDefinition() )
        leadingStrangeParticleHasChanged = false;
      if( leadingStrangeParticleHasChanged &&
          ( targetParticle.GetDefinition() == leadingStrangeParticle.GetDefinition() ) )
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
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
        }
        else
        {
          currentParticle.SetDefinitionAndUpdateE( leadingStrangeParticle.GetDefinition() );
          incidentHasChanged = false;
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
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
      G4int(targetNucleus.GetZ()) - protonsInFinalState);
    G4int NinNucleus = std::max(0,
      G4int(targetNucleus.GetN()-targetNucleus.GetZ()) - neutronsInFinalState);

    pseudoParticle[3].SetMomentum( 0.0, 0.0, pOriginal*GeV );
    pseudoParticle[3].SetMass( mOriginal*GeV );
    pseudoParticle[3].SetTotalEnergy(
     std::sqrt( pOriginal*pOriginal + mOriginal*mOriginal )*GeV );
    
    G4ParticleDefinition * aOrgDef = modifiedOriginal.GetDefinition();
    G4int diff = 0;
    if(aOrgDef == G4Proton::Proton() || aOrgDef == G4Neutron::Neutron() )  diff = 1;
    if(numberofFinalStateNucleons == 1) diff = 0;
    pseudoParticle[4].SetMomentum( 0.0, 0.0, 0.0 );
    pseudoParticle[4].SetMass( protonMass*(numberofFinalStateNucleons-diff)*MeV );
    pseudoParticle[4].SetTotalEnergy( protonMass*(numberofFinalStateNucleons-diff)*MeV );
    
    G4double theoreticalKinetic =
      pseudoParticle[3].GetTotalEnergy()/MeV +
      pseudoParticle[4].GetTotalEnergy()/MeV -
      currentParticle.GetMass()/MeV -
      targetParticle.GetMass()/MeV;
    
    G4double simulatedKinetic =
      currentParticle.GetKineticEnergy()/MeV + targetParticle.GetKineticEnergy()/MeV;
    
    pseudoParticle[5] = pseudoParticle[3] + pseudoParticle[4];
    pseudoParticle[3].Lorentz( pseudoParticle[3], pseudoParticle[5] );
    pseudoParticle[4].Lorentz( pseudoParticle[4], pseudoParticle[5] );
    
    pseudoParticle[7].SetZero();
    pseudoParticle[7] = pseudoParticle[7] + currentParticle;
    pseudoParticle[7] = pseudoParticle[7] + targetParticle;

    for( i=0; i<vecLen; ++i )
    {
      pseudoParticle[7] = pseudoParticle[7] + *vec[i];
      simulatedKinetic += vec[i]->GetKineticEnergy()/MeV;
      theoreticalKinetic -= vec[i]->GetMass()/MeV;
    }

    if( vecLen <= 16 && vecLen > 0 )
    {
      // must create a new set of ReactionProducts here because GenerateNBody will
      // modify the momenta for the particles, and we don't want to do this
      //
      G4ReactionProduct tempR[130];
      tempR[0] = currentParticle;
      tempR[1] = targetParticle;
      for( i=0; i<vecLen; ++i )tempR[i+2] = *vec[i];
      G4FastVector<G4ReactionProduct,256> tempV;
      tempV.Initialize( vecLen+2 );
      G4int tempLen = 0;
      for( i=0; i<vecLen+2; ++i )tempV.SetElement( tempLen++, &tempR[i] );
      constantCrossSection = true;

      wgt = GenerateNBodyEvent( pseudoParticle[3].GetTotalEnergy()/MeV+
                                pseudoParticle[4].GetTotalEnergy()/MeV,
                                constantCrossSection, tempV, tempLen );
      if (wgt == -1) {
        G4double Qvalue = 0;
        for (i = 0; i < tempLen; i++) Qvalue += tempV[i]->GetMass();
        wgt = GenerateNBodyEvent( Qvalue/MeV,
                                  constantCrossSection, tempV, tempLen );
      }
      if(wgt>-.5)
      {
        theoreticalKinetic = 0.0;
        for( i=0; i<tempLen; ++i )
        {
          pseudoParticle[6].Lorentz( *tempV[i], pseudoParticle[4] );
          theoreticalKinetic += pseudoParticle[6].GetKineticEnergy()/MeV;
        }
      }
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    }
    //
    //  Make sure, that the kinetic energies are correct
    //
    if( simulatedKinetic != 0.0 )
    {
      wgt = (theoreticalKinetic)/simulatedKinetic;
      theoreticalKinetic = currentParticle.GetKineticEnergy()/MeV * wgt;
      simulatedKinetic = theoreticalKinetic;
      currentParticle.SetKineticEnergy( theoreticalKinetic*MeV );
      pp = currentParticle.GetTotalMomentum()/MeV;
      pp1 = currentParticle.GetMomentum().mag()/MeV;
      if( pp1 < 1.0e-6*GeV ) {
        G4ThreeVector iso = Isotropic(pp);
        currentParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
      }
      theoreticalKinetic = targetParticle.GetKineticEnergy()/MeV * wgt;
      targetParticle.SetKineticEnergy( theoreticalKinetic*MeV );
      simulatedKinetic += theoreticalKinetic;
      pp = targetParticle.GetTotalMomentum()/MeV;
      pp1 = targetParticle.GetMomentum().mag()/MeV;
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
      if( pp1 < 1.0e-6*GeV ) {
        G4ThreeVector iso = Isotropic(pp);
        targetParticle.SetMomentum( iso.x(), iso.y(), iso.z() );
      } else {
        targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
      }

      for( i=0; i<vecLen; ++i ) {
        theoreticalKinetic = vec[i]->GetKineticEnergy()/MeV * wgt;
        simulatedKinetic += theoreticalKinetic;
        vec[i]->SetKineticEnergy( theoreticalKinetic*MeV );
        pp = vec[i]->GetTotalMomentum()/MeV;
        pp1 = vec[i]->GetMomentum().mag()/MeV;
        if( pp1 < 1.0e-6*GeV ) {
          G4ThreeVector iso = Isotropic(pp);
          vec[i]->SetMomentum( iso.x(), iso.y(), iso.z() );
        } else {
          vec[i]->SetMomentum( vec[i]->GetMomentum() * (pp/pp1) );
	}
      }
    }
    // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);

    Rotate( numberofFinalStateNucleons, pseudoParticle[3].GetMomentum(),
            modifiedOriginal, originalIncident, targetNucleus,
            currentParticle, targetParticle, vec, vecLen );
      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    //
    // add black track particles
    // the total number of particles produced is restricted to 198
    // this may have influence on very high energies
    //
    if( atomicWeight >= 1.5 )
    {
      // npnb is number of proton/neutron black track particles
      // ndta is the number of deuterons, tritons, and alphas produced
      // epnb is the kinetic energy available for proton/neutron black track particles
      // edta is the kinetic energy available for deuteron/triton/alpha particles
      //
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

      const G4double pnCutOff = 0.001;
      const G4double dtaCutOff = 0.001;
      const G4double kineticMinimum = 1.e-6;
      const G4double kineticFactor = -0.010;
      G4double sprob = 0.0; // sprob = probability of self-absorption in heavy molecules
      const G4double ekIncident = originalIncident->GetKineticEnergy()/GeV;
      if( ekIncident >= 5.0 )sprob = std::min( 1.0, 0.6*std::log(ekIncident-4.0) );
      if( epnb >= pnCutOff )
      {
        npnb = G4Poisson((1.5+1.25*numberofFinalStateNucleons)*epnb/(epnb+edta));
        if( numberofFinalStateNucleons + npnb > atomicWeight )
          npnb = G4int(atomicWeight+0.00001 - numberofFinalStateNucleons);
        npnb = std::min( npnb, 127-vecLen );
      }
      if( edta >= dtaCutOff )
      {
        ndta = G4Poisson( (1.5+1.25*numberofFinalStateNucleons)*edta/(epnb+edta) );
        ndta = std::min( ndta, 127-vecLen );
      }
      if (npnb == 0 && ndta == 0) npnb = 1;

      // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);

      AddBlackTrackParticles(epnb, npnb, edta, ndta, sprob, kineticMinimum, 
                             kineticFactor, modifiedOriginal,
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
         1.0-500.0*std::exp(-ekOriginal/0.04)*std::log(G4UniformRand()) );
  else
    currentParticle.SetTOF( 1.0 );
  return true;

}
 
 /* end of file */
