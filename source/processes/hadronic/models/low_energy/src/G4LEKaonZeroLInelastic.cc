// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEKaonZeroLInelastic.cc,v 1.1 1999-01-07 16:12:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy KaonZeroLong Inelastic Process
 // J.L. Chuma, TRIUMF, 11-Feb-1997
 // Last modified: 27-Mar-1997
 // Modified by J.L.Chuma 30-Apr-97: added originalTarget for CalculateMomenta
 
#include "G4LEKaonZeroLInelastic.hh"
#include "Randomize.hh"
 
 G4VParticleChange *
  G4LEKaonZeroLInelastic::ApplyYourself( const G4Track &aTrack,
                                         G4Nucleus &targetNucleus )
  {
    theParticleChange.Initialize( aTrack );
    
    const G4DynamicParticle *originalIncident = aTrack.GetDynamicParticle();
    //
    // create the target particle
    //
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    
    if( verboseLevel > 1 )
    {
      G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEKaonZeroLInelastic::ApplyYourself called" << endl;    
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
      G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << endl;
    }
    //
    // Fermi motion and evaporation
    // As of Geant3, the Fermi energy calculation had not been Done
    //
    G4double ek = originalIncident->GetKineticEnergy()/MeV;
    G4double amas = originalIncident->GetDefinition()->GetPDGMass()/MeV;
    G4ReactionProduct modifiedOriginal;
    modifiedOriginal = *originalIncident;
    
    G4double tkin = targetNucleus.Cinema( ek );
    ek += tkin;
    modifiedOriginal.SetKineticEnergy( ek*MeV );
    G4double et = ek + amas;
    G4double p = sqrt( abs((et-amas)*(et+amas)) );
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
    modifiedOriginal.SetKineticEnergy( ek*MeV );
    et = ek + amas;
    p = sqrt( abs((et-amas)*(et+amas)) );
    pp = modifiedOriginal.GetMomentum().mag()/MeV;
    if( pp > 0.0 )
    {
      G4ThreeVector momentum = modifiedOriginal.GetMomentum();
      modifiedOriginal.SetMomentum( momentum * (p/pp) );
    }
    G4ReactionProduct currentParticle = modifiedOriginal;
    G4ReactionProduct targetParticle;
    targetParticle = *originalTarget;
    currentParticle.SetSide( 1 ); // incident always goes in forward hemisphere
    targetParticle.SetSide( -1 );  // target always goes in backward hemisphere
    G4bool incidentHasChanged = false;
    G4bool targetHasChanged = false;
    G4bool quasiElastic = false;
    G4FastVector<G4ReactionProduct,128> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    const G4double cutOff = 0.1;
    if( currentParticle.GetKineticEnergy()/MeV > cutOff )
      Cascade( vec, vecLen,
               originalIncident, currentParticle, targetParticle,
               incidentHasChanged, targetHasChanged, quasiElastic );
    
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
 
 void
  G4LEKaonZeroLInelastic::Cascade(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int& vecLen,
   const G4DynamicParticle *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASK0B by H. Fesefeldt (13-Sep-1987)
    //
    // K0Long undergoes interaction with nucleon within a nucleus.  Check if it is
    // energetically possible to produce pions/kaons.  In not, assume nuclear excitation
    // occurs and input particle is degraded in energy. No other particles are produced.
    // If reaction is possible, find the correct number of pions/protons/neutrons
    // produced using an interpolation to multiplicity data.  Replace some pions or
    // protons/neutrons by kaons or strange baryons according to the average
    // multiplicity per Inelastic reaction.
    //
    const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass()/MeV;
    const G4double etOriginal = originalIncident->GetTotalEnergy()/MeV;
    const G4double pOriginal = originalIncident->GetTotalMomentum()/MeV;
    const G4double targetMass = targetParticle.GetMass()/MeV;
    G4double centerofmassEnergy = sqrt( mOriginal*mOriginal +
                                        targetMass*targetMass +
                                        2.0*targetMass*etOriginal );
    G4double availableEnergy = centerofmassEnergy-(targetMass+mOriginal);

    static G4bool first = true;
    const G4int numMul = 1200;
    const G4int numSec = 60;
    static G4double protmul[numMul], protnorm[numSec]; // proton constants
    static G4double neutmul[numMul], neutnorm[numSec]; // neutron constants
    // np = number of pi+, nm = number of pi-, nz = number of pi0
    G4int counter, nt=0, np=0, nm=0, nz=0;
    const G4double c = 1.25;    
    const G4double b[] = { 0.7, 0.7 };
    if( first )       // compute normalization constants, this will only be Done once
    {
      first = false;
      G4int i;
      for( i=0; i<numMul; ++i )protmul[i] = 0.0;
      for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
      counter = -1;
      for( np=0; np<numSec/3; ++np )
      {
        for( nm=max(0,np-2); nm<=np; ++nm )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec )
              {
                protmul[counter] = Pmltpc(np,nm,nz,nt,b[0],c);
                protnorm[nt-1] += protmul[counter];
              }
            }
          }
        }
      }
      for( i=0; i<numMul; ++i )neutmul[i] = 0.0;
      for( i=0; i<numSec; ++i )neutnorm[i] = 0.0;
      counter = -1;
      for( np=0; np<(numSec/3); ++np )
      {
        for( nm=max(0,np-1); nm<=(np+1); ++nm )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec )
              {
                neutmul[counter] = Pmltpc(np,nm,nz,nt,b[1],c);
                neutnorm[nt-1] += neutmul[counter];
              }
            }
          }
        }
      }
      for( i=0; i<numSec; ++i )
      {
        if( protnorm[i] > 0.0 )protnorm[i] = 1.0/protnorm[i];
        if( neutnorm[i] > 0.0 )neutnorm[i] = 1.0/neutnorm[i];
      }
    }   // end of initialization
    
    const G4double expxu = 82.;           // upper bound for arg. of exp
    const G4double expxl = -expxu;        // lower bound for arg. of exp
    G4ParticleDefinition *aKaonPlus = G4KaonPlus::KaonPlus();
    G4ParticleDefinition *aKaonMinus = G4KaonMinus::KaonMinus();
    G4ParticleDefinition *aKaonZS = G4KaonZeroShort::KaonZeroShort();
    G4ParticleDefinition *aKaonZL = G4KaonZeroLong::KaonZeroLong();
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aPiPlus = G4PionPlus::PionPlus();
    G4ParticleDefinition *aPiMinus = G4PionMinus::PionMinus();
    G4ParticleDefinition *aPiZero = G4PionZero::PionZero();
    G4ParticleDefinition *aLambda = G4Lambda::Lambda();
    G4ParticleDefinition *aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition *aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition *aSigmaZero = G4SigmaZero::SigmaZero();
    const G4double cech[] = {1.,1.,1.,0.70,0.60,0.55,0.35,0.25,0.18,0.15};
    G4int iplab = min( 9.0, 5.0*pOriginal*MeV/GeV );
    if( (pOriginal*MeV/GeV <= 2.0) && (G4UniformRand() < cech[iplab]) )
    {
      np = nm = nz = nt = 0;
      iplab = min( 19.0, pOriginal*MeV/GeV*10.0 );
      const G4double cnk0[] = {0.17,0.18,0.17,0.24,0.26,0.20,0.22,0.21,0.34,0.45,
                               0.58,0.55,0.36,0.29,0.29,0.32,0.32,0.33,0.33,0.33};
      if( G4UniformRand() > cnk0[iplab] )
      {
        G4double ran = G4UniformRand();
        if( ran < 0.25 )         // k0Long n --> pi- s+
        {
          if( targetParticle.GetDefinition() == aNeutron )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
        }
        else if( ran < 0.50 )  // k0Long p --> pi+ s0  or  k0Long n --> pi0 s0
        {
          if( targetParticle.GetDefinition() == aNeutron )
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
          else
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
          targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
        else if( ran < 0.75 )  // k0Long n --> pi+ s-
        {
          if( targetParticle.GetDefinition() == aNeutron )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
        }
        else                   // k0Long p --> pi+ L  or  k0Long n --> pi0 L
        {
          if( targetParticle.GetDefinition() == aNeutron )
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
          else
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
          targetParticle.SetDefinitionAndUpdateE( aLambda );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
      }
      else   // ran <= cnk0
      {
        quasiElastic = true;
        if( targetParticle.GetDefinition() == aNeutron )
        {
          currentParticle.SetDefinitionAndUpdateE( aKaonMinus );
          targetParticle.SetDefinitionAndUpdateE( aProton );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
      }
    }
    else  // (pOriginal > 2.0*GeV) || (random number >= cech[iplab])
    {
      if( availableEnergy < aPiPlus->GetPDGMass()/MeV )
      {
        quasiElastic = true;
        return;
      }
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, test, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        counter = -1;
        for( np=0; (np<numSec/3) && (ran>=excs); ++np )
        {
          for( nm=max(0,np-2); nm<=np && ran>=excs; ++nm )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( nt>0 && nt<=numSec )
                {
                  test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
                  if( fabs(dum) < 1.0 )
                  {
                    if( test >= 1.0e-10 )excs += dum*test;
                  }
                  else
                    excs += dum*test;
                }
              }
            }
          }
        }
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nm--; nz--;
        switch( np-nm )
        {
         case 1:
           if( G4UniformRand() < 0.5 )
           {
             currentParticle.SetDefinitionAndUpdateE( aKaonMinus );
             incidentHasChanged = true;
           }
           else
           {
             targetParticle.SetDefinitionAndUpdateE( aNeutron );
             targetHasChanged = true;
           }
         case 0:
           break;
         default:
           currentParticle.SetDefinitionAndUpdateE( aKaonMinus );
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           incidentHasChanged = true;
           targetHasChanged = true;
           break;
        }
      }
      else  // target must be a neutron
      {
        counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nm=max(0,np-1); nm<=(np+1) && ran>=excs; ++nm )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( nt>0 && nt<=numSec )
                {
                  test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
                  if( fabs(dum) < 1.0 )
                  {
                    if( test >= 1.0e-10 )excs += dum*test;
                  }
                  else
                    excs += dum*test;
                }
              }
            }
          }
        }
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nm--; nz--;
        switch( np-nm )
        {
         case 0:
           currentParticle.SetDefinitionAndUpdateE( aKaonMinus );
           targetParticle.SetDefinitionAndUpdateE( aProton );
           incidentHasChanged = true;
           targetHasChanged = true;
           break;
         case 1:
           currentParticle.SetDefinitionAndUpdateE( aKaonMinus );
           incidentHasChanged = true;
           break;
         default:
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
           break;
        }
      }
      if( G4UniformRand() >= 0.5 )
      {
        if( currentParticle.GetDefinition() == aKaonMinus &&
            targetParticle.GetDefinition() == aNeutron )
        {
          ran = G4UniformRand();
          if( ran < 0.68 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aLambda );
          }
          else if( ran < 0.84 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
          }
          else
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
          }
        }
        else if( (currentParticle.GetDefinition() == aKaonZS ||
                  currentParticle.GetDefinition() == aKaonZL ) &&
                 targetParticle.GetDefinition() == aProton )
        {
          ran = G4UniformRand();
          if( ran < 0.68 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
            targetParticle.SetDefinitionAndUpdateE( aLambda );
          }
          else if( ran < 0.84 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
          }
          else
          {
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
          }
        }
        else
        {
          ran = G4UniformRand();
          if( ran < 0.67 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aLambda );
          }
          else if( ran < 0.78 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
          }
          else if( ran < 0.89 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
          }
          else
          {
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
          }
        }
        incidentHasChanged = true;
        targetHasChanged = true;
      }
    }
    if( currentParticle.GetDefinition() == aKaonZL )
    {
      if( G4UniformRand() >= 0.5 )
      {
        currentParticle.SetDefinitionAndUpdateE( aKaonZS );
        incidentHasChanged = true;
      }
    }
    if( targetParticle.GetDefinition() == aKaonZL )
    {
      if( G4UniformRand() >= 0.5 )
      {
        targetParticle.SetDefinitionAndUpdateE( aKaonZS );
        targetHasChanged = true;
      }
    }
    SetUpPions( np, nm, nz, vec, vecLen );
    return;
  }

 /* end of file */
 
