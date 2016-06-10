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
// $Id: G4RPGStrangeProduction.cc 94214 2015-11-09 08:18:05Z gcosmo $
//
 
#include <iostream>
#include <signal.h>

#include "G4RPGStrangeProduction.hh"
#include "G4Log.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadReentrentException.hh"

G4RPGStrangeProduction::G4RPGStrangeProduction()
  : G4RPGReaction() {}


G4bool G4RPGStrangeProduction::
ReactionStage(const G4HadProjectile* /*originalIncident*/,
              G4ReactionProduct& modifiedOriginal,
              G4bool& incidentHasChanged,
              const G4DynamicParticle* originalTarget,
              G4ReactionProduct& targetParticle,
              G4bool& targetHasChanged,
              const G4Nucleus& /*targetNucleus*/,
              G4ReactionProduct& currentParticle,
              G4FastVector<G4ReactionProduct,256>& vec,
              G4int& vecLen,
              G4bool /*leadFlag*/,
              G4ReactionProduct& /*leadingStrangeParticle*/)
{
  // Derived from H. Fesefeldt's original FORTRAN code STPAIR
  //
  // Choose charge combinations K+ K-, K+ K0B, K0 K0B, K0 K-,
  //                            K+ Y0, K0 Y+,  K0 Y-
  // For antibaryon induced reactions half of the cross sections KB YB
  // pairs are produced.  Charge is not conserved, no experimental data available
  // for exclusive reactions, therefore some average behaviour assumed.
  // The ratio L/SIGMA is taken as 3:1 (from experimental low energy)
  //

  if( vecLen == 0 )return true;
  //
  // the following protects against annihilation processes
  //
  if( currentParticle.GetMass() == 0.0 || targetParticle.GetMass() == 0.0 )return true;
    
  const G4double etOriginal = modifiedOriginal.GetTotalEnergy()/GeV;
  const G4double mOriginal = modifiedOriginal.GetDefinition()->GetPDGMass()/GeV;
  G4double targetMass = originalTarget->GetDefinition()->GetPDGMass()/GeV;
  G4double centerofmassEnergy = std::sqrt( mOriginal*mOriginal +
                                      targetMass*targetMass +
                                      2.0*targetMass*etOriginal );  // GeV
  G4double currentMass = currentParticle.GetMass()/GeV;
  G4double availableEnergy = centerofmassEnergy-(targetMass+currentMass);
  if( availableEnergy <= 1.0 )return true;
    
  G4ParticleDefinition *aProton = G4Proton::Proton();
  G4ParticleDefinition *anAntiProton = G4AntiProton::AntiProton();
  G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
  G4ParticleDefinition *anAntiNeutron = G4AntiNeutron::AntiNeutron();
  G4ParticleDefinition *aSigmaMinus = G4SigmaMinus::SigmaMinus();
  G4ParticleDefinition *aSigmaPlus = G4SigmaPlus::SigmaPlus();
  G4ParticleDefinition *aSigmaZero = G4SigmaZero::SigmaZero();
  G4ParticleDefinition *anAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinus();
  G4ParticleDefinition *anAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlus();
  G4ParticleDefinition *anAntiSigmaZero = G4AntiSigmaZero::AntiSigmaZero();
  G4ParticleDefinition *aKaonMinus = G4KaonMinus::KaonMinus();
  G4ParticleDefinition *aKaonPlus = G4KaonPlus::KaonPlus();
  G4ParticleDefinition *aKaonZL = G4KaonZeroLong::KaonZeroLong();
  G4ParticleDefinition *aKaonZS = G4KaonZeroShort::KaonZeroShort();
  G4ParticleDefinition *aLambda = G4Lambda::Lambda();
  G4ParticleDefinition *anAntiLambda = G4AntiLambda::AntiLambda();
    
  const G4double protonMass = aProton->GetPDGMass()/GeV;
  const G4double sigmaMinusMass = aSigmaMinus->GetPDGMass()/GeV;
  //
  // determine the center of mass energy bin
  //
  const G4double avrs[] = {3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.};

  G4int ibin, i3, i4;
  G4double avk, avy, avn, ran;
  G4int i = 1;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
  while ((i<12) && (centerofmassEnergy>avrs[i]) ) {   /* Loop checking, 01.09.2015, D.Wright */
    i++;
    loop++;
    if (loop > 1000) {
      G4Exception("G4RPGStrangeProduction::ReactionStage()", "HAD_RPG_100", JustWarning, ed);
      break;
    }
  }


  if( i == 12 )
    ibin = 11;
  else
    ibin = i;
  //
  // the fortran code chooses a random replacement of produced kaons
  //  but does not take into account charge conservation 
  //
  if (vecLen == 1) { // we know that vecLen > 0
    i3 = 0;
    i4 = 1;   // note that we will be adding a new secondary particle in this case only
  } else {    // otherwise  0 <= i3,i4 < vecLen
    G4double rnd = G4UniformRand();
    while (rnd == 1.0) rnd = G4UniformRand();  /* Loop checking, 01.09.2015, D.Wright */
    i4 = i3 = G4int(vecLen*rnd);

    while(i3 == i4) {  /* Loop checking, 01.09.2015, D.Wright */
      rnd = G4UniformRand();
      while( rnd == 1.0 ) rnd = G4UniformRand();  /* Loop checking, 01.09.2015, D.Wright */
      i4 = G4int(vecLen*rnd);
    }
  }

  // use linear interpolation or extrapolation by y=centerofmassEnergy*x+b
  //
  const G4double avkkb[] = { 0.0015, 0.005, 0.012, 0.0285, 0.0525, 0.075,
                             0.0975, 0.123, 0.28,  0.398,  0.495,  0.573 };
  const G4double avky[] = { 0.005, 0.03,  0.064, 0.095, 0.115, 0.13,
                            0.145, 0.155, 0.20,  0.205, 0.210, 0.212 };
  const G4double avnnb[] = { 0.00001, 0.0001, 0.0006, 0.0025, 0.01, 0.02,
                             0.04,    0.05,   0.12,   0.15,   0.18, 0.20 };
    
  avk = (G4Log(avkkb[ibin])-G4Log(avkkb[ibin-1]))*(centerofmassEnergy-avrs[ibin-1])
    /(avrs[ibin]-avrs[ibin-1]) + G4Log(avkkb[ibin-1]);
  avk = G4Exp(avk);
    
  avy = (G4Log(avky[ibin])-G4Log(avky[ibin-1]))*(centerofmassEnergy-avrs[ibin-1])
    /(avrs[ibin]-avrs[ibin-1]) + G4Log(avky[ibin-1]);
  avy = G4Exp(avy);
    
  avn = (G4Log(avnnb[ibin])-G4Log(avnnb[ibin-1]))*(centerofmassEnergy-avrs[ibin-1])
    /(avrs[ibin]-avrs[ibin-1]) + G4Log(avnnb[ibin-1]);
  avn = G4Exp(avn);
    
  if( avk+avy+avn <= 0.0 )return true;
    
  if( currentMass < protonMass )avy /= 2.0;
  if( targetMass < protonMass )avy = 0.0;
  avy += avk+avn;
  avk += avn;
  ran = G4UniformRand();
  if(  ran < avn )
  {
    if( availableEnergy < 2.0 )return true;
    if( vecLen == 1 )                              // add a new secondary
    {
      G4ReactionProduct *p1 = new G4ReactionProduct;
      if( G4UniformRand() < 0.5 )
      {
        vec[0]->SetDefinition( aNeutron );
        p1->SetDefinition( anAntiNeutron );
        (G4UniformRand() < 0.5) ? p1->SetSide( -1 ) : p1->SetSide( 1 );
        vec[0]->SetMayBeKilled(false);
        p1->SetMayBeKilled(false);
      }
      else
      {
        vec[0]->SetDefinition( aProton );
        p1->SetDefinition( anAntiProton );
        (G4UniformRand() < 0.5) ? p1->SetSide( -1 ) : p1->SetSide( 1 );
        vec[0]->SetMayBeKilled(false);
        p1->SetMayBeKilled(false);
      }
      vec.SetElement( vecLen++, p1 );
    // DEBUGGING --> DumpFrames::DumpFrame(vec, vecLen);
    }
    else
    {                                             // replace two secondaries
      if( G4UniformRand() < 0.5 )
      {
        vec[i3]->SetDefinition( aNeutron );
        vec[i4]->SetDefinition( anAntiNeutron );
        vec[i3]->SetMayBeKilled(false);
        vec[i4]->SetMayBeKilled(false);
      }
      else
      {
        vec[i3]->SetDefinition( aProton );
        vec[i4]->SetDefinition( anAntiProton );
        vec[i3]->SetMayBeKilled(false);
        vec[i4]->SetMayBeKilled(false);
      }
    }
  }
  else if( ran < avk )
  {
    if( availableEnergy < 1.0 )return true;
      
    const G4double kkb[] = { 0.2500, 0.3750, 0.5000, 0.5625, 0.6250,
                               0.6875, 0.7500, 0.8750, 1.000 };
    const G4int ipakkb1[] = { 10, 10, 10, 11, 11, 12, 12, 11, 12 };
    const G4int ipakkb2[] = { 13, 11, 12, 11, 12, 11, 12, 13, 13 };
    ran = G4UniformRand();
    i = 0;

    loop = 0;
    G4ExceptionDescription eda;
    eda << " While count exceeded " << G4endl;
    while( (i<9) && (ran>=kkb[i]) ) {  /* Loop checking, 01.09.2015, D.Wright */
      ++i;
      loop++;
      if (loop > 1000) {
        G4Exception("G4RPGStrangeProduction::ReactionStage()", "HAD_RPG_100", JustWarning, eda);
        break;
      }
    }

    if( i == 9 )return true;
    //
    // ipakkb[] = { 10,13, 10,11, 10,12, 11,11, 11,12, 12,11, 12,12, 11,13, 12,13 };
    // charge       +  -   +  0   +  0   0  0   0  0   0  0   0  0   0  -   0  -
    //
    switch( ipakkb1[i] )
    {
     case 10:
       vec[i3]->SetDefinition( aKaonPlus );
       vec[i3]->SetMayBeKilled(false);
       break;
     case 11:
       vec[i3]->SetDefinition( aKaonZS );
       vec[i3]->SetMayBeKilled(false);
       break;
     case 12:
       vec[i3]->SetDefinition( aKaonZL );
       vec[i3]->SetMayBeKilled(false);
       break;
    }

    if( vecLen == 1 )                          // add a secondary
    {
      G4ReactionProduct *p1 = new G4ReactionProduct;
      switch( ipakkb2[i] )
      {
       case 11:
         p1->SetDefinition( aKaonZS );
         p1->SetMayBeKilled(false);
         break;
       case 12:
         p1->SetDefinition( aKaonZL );
         p1->SetMayBeKilled(false);
         break;
       case 13:
         p1->SetDefinition( aKaonMinus );
	 p1->SetMayBeKilled(false);
         break;
      }
      (G4UniformRand() < 0.5) ? p1->SetSide( -1 ) : p1->SetSide( 1 );
      vec.SetElement( vecLen++, p1 );

    }
    else                                        // replace
    {
      switch( ipakkb2[i] )
      {
       case 11:
         vec[i4]->SetDefinition( aKaonZS );
         vec[i4]->SetMayBeKilled(false);
         break;
       case 12:
         vec[i4]->SetDefinition( aKaonZL );
         vec[i4]->SetMayBeKilled(false);
         break;
       case 13:
         vec[i4]->SetDefinition( aKaonMinus );
	 vec[i4]->SetMayBeKilled(false);
         break;
      }
    }

  } else if( ran < avy ) {
    if( availableEnergy < 1.6 )return true;
      
    const G4double ky[] = { 0.200, 0.300, 0.400, 0.550, 0.625, 0.700,
                            0.800, 0.850, 0.900, 0.950, 0.975, 1.000 };
    const G4int ipaky1[] = { 18, 18, 18, 20, 20, 20, 21, 21, 21, 22, 22, 22 };
    const G4int ipaky2[] = { 10, 11, 12, 10, 11, 12, 10, 11, 12, 10, 11, 12 };
    const G4int ipakyb1[] = { 19, 19, 19, 23, 23, 23, 24, 24, 24, 25, 25, 25 };
    const G4int ipakyb2[] = { 13, 12, 11, 13, 12, 11, 13, 12, 11, 13, 12, 11 };
    ran = G4UniformRand();
    i = 0;

    loop = 0;
    G4ExceptionDescription edb;
    edb << " While count exceeded " << G4endl;
    while( (i<12) && (ran>ky[i]) ) {  /* Loop checking, 01.09.2015, D.Wright */
      ++i;
      loop++;
      if (loop > 1000) {
        G4Exception("G4RPGStrangeProduction::ReactionStage()", "HAD_RPG_100", JustWarning, edb);
        break;
      }
    }

    if( i == 12 )return true;
    if ( (currentMass<protonMass) || (G4UniformRand()<0.5) ) {
      // ipaky[] = { 18,10, 18,11, 18,12, 20,10, 20,11, 20,12,
      //             0  +   0  0   0  0   +  +   +  0   +  0
      //
      //             21,10, 21,11, 21,12, 22,10, 22,11, 22,12 }
      //             0  +   0  0   0  0   -  +   -  0   -  0
      switch( ipaky1[i] )
      {
       case 18:
         targetParticle.SetDefinition( aLambda );
         break;
       case 20:
         targetParticle.SetDefinition( aSigmaPlus );
         break;
       case 21:
         targetParticle.SetDefinition( aSigmaZero );
         break;
       case 22:
         targetParticle.SetDefinition( aSigmaMinus );
         break;
      }
      targetHasChanged = true;
      switch( ipaky2[i] )
      {
       case 10:
         vec[i3]->SetDefinition( aKaonPlus ); 
         vec[i3]->SetMayBeKilled(false);
         break;
       case 11:
         vec[i3]->SetDefinition( aKaonZS );
         vec[i3]->SetMayBeKilled(false);
         break;
       case 12:
         vec[i3]->SetDefinition( aKaonZL );
         vec[i3]->SetMayBeKilled(false);
         break;
      }

    } else {  // (currentMass >= protonMass) && (G4UniformRand() >= 0.5)
      // ipakyb[] = { 19,13, 19,12, 19,11, 23,13, 23,12, 23,11,
      //              24,13, 24,12, 24,11, 25,13, 25,12, 25,11 };
      if ( (currentParticle.GetDefinition() == anAntiProton) ||
           (currentParticle.GetDefinition() == anAntiNeutron) ||
           (currentParticle.GetDefinition() == anAntiLambda) ||
           (currentMass > sigmaMinusMass) ) {
        switch( ipakyb1[i] )
        {
         case 19:
           currentParticle.SetDefinitionAndUpdateE( anAntiLambda );
           break;
         case 23:
           currentParticle.SetDefinitionAndUpdateE( anAntiSigmaPlus );
           break;
         case 24:
           currentParticle.SetDefinitionAndUpdateE( anAntiSigmaZero );
           break;
         case 25:
           currentParticle.SetDefinitionAndUpdateE( anAntiSigmaMinus );
           break;
        }
        incidentHasChanged = true;
        switch( ipakyb2[i] )
        {
         case 11:
           vec[i3]->SetDefinition( aKaonZS ); 
           vec[i3]->SetMayBeKilled(false);
           break;
         case 12:
           vec[i3]->SetDefinition( aKaonZL );
           vec[i3]->SetMayBeKilled(false);
           break;
         case 13:
           vec[i3]->SetDefinition( aKaonMinus );
           vec[i3]->SetMayBeKilled(false);
           break;
        }

      } else {
        switch( ipaky1[i] )
        {
         case 18:
           currentParticle.SetDefinitionAndUpdateE( aLambda );
           break;
         case 20:
           currentParticle.SetDefinitionAndUpdateE( aSigmaPlus );
           break;
         case 21:
           currentParticle.SetDefinitionAndUpdateE( aSigmaZero );
           break;
         case 22:
           currentParticle.SetDefinitionAndUpdateE( aSigmaMinus );
           break;
        }
        incidentHasChanged = true;
        switch( ipaky2[i] )
        {
         case 10:
           vec[i3]->SetDefinition( aKaonPlus ); 
           vec[i3]->SetMayBeKilled(false);
           break;
         case 11:
           vec[i3]->SetDefinition( aKaonZS );
           vec[i3]->SetMayBeKilled(false);
           break;
         case 12:
           vec[i3]->SetDefinition( aKaonZL );
           vec[i3]->SetMayBeKilled(false);
           break;
        }
      }
    }
  }
  else return true;

  //
  //  check the available energy
  //   if there is not enough energy for kkb/ky pair production
  //   then reduce the number of secondary particles 
  //  NOTE:
  //        the number of secondaries may have been changed
  //        the incident and/or target particles may have changed
  //        charge conservation is ignored (as well as strangness conservation)
  //
  currentMass = currentParticle.GetMass()/GeV;
  targetMass = targetParticle.GetMass()/GeV;
    
  G4double energyCheck = centerofmassEnergy-(currentMass+targetMass);
  for( i=0; i<vecLen; ++i )
  {
    energyCheck -= vec[i]->GetMass()/GeV;
    if( energyCheck < 0.0 )      // chop off the secondary List
    {
      vecLen = std::max( 0, --i ); // looks like a memory leak @@@@@@@@@@@@
      G4int j;
      for(j=i; j<vecLen; j++) delete vec[j];
      break;
    }
  }

  return true;
}

 
 /* end of file */
