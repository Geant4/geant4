// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEKaonMinusInelastic.cc,v 1.2 1999-12-15 14:53:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy KaonMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 12-Feb-1997
 // Last modified: 27-Mar-1997
 // J.P.Wellisch 23-Apr-97: bug-hunting (missing initialization of np,nm,nz fixed)
 // Modified by J.L.Chuma 30-Apr-97: added originalTarget for CalculateMomenta
 
#include "G4LEKaonMinusInelastic.hh"
 #include "Randomize.hh"

 G4VParticleChange *
  G4LEKaonMinusInelastic::ApplyYourself( const G4Track &aTrack,
                                         G4Nucleus &targetNucleus )
  {
    theParticleChange.Initialize( aTrack );
    
    const G4DynamicParticle *originalIncident = aTrack.GetDynamicParticle();
    if (originalIncident->GetKineticEnergy()<= 0.1*MeV) return &theParticleChange;

    // create the target particle
    
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    G4double targetMass = originalTarget->GetDefinition()->GetPDGMass();
    G4ReactionProduct targetParticle( originalTarget->GetDefinition() );
    
    if( verboseLevel > 1 )
    {
      G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEKaonMinusInelastic::ApplyYourself called" << G4endl;
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy() << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
      G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << G4endl;
    }
    G4ReactionProduct currentParticle( originalIncident->GetDefinition() );
    currentParticle.SetMomentum( originalIncident->GetMomentum() );
    currentParticle.SetKineticEnergy( originalIncident->GetKineticEnergy() );
    
    // Fermi motion and evaporation
    // As of Geant3, the Fermi energy calculation had not been Done
    
    G4double ek = originalIncident->GetKineticEnergy();
    G4double amas = originalIncident->GetDefinition()->GetPDGMass();
    
    G4double tkin = targetNucleus.Cinema( ek );
    ek += tkin;
    currentParticle.SetKineticEnergy( ek );
    G4double et = ek + amas;
    G4double p = sqrt( abs((et-amas)*(et+amas)) );
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
    p = sqrt( abs((et-amas)*(et+amas)) );
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
    G4FastVector<G4ReactionProduct,128> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    const G4double cutOff = 0.1*MeV;
    if( currentParticle.GetKineticEnergy() > cutOff )
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
  G4LEKaonMinusInelastic::Cascade(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int& vecLen,
   const G4DynamicParticle *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASKM by H. Fesefeldt (13-Sep-1987)
    //
    // K- undergoes interaction with nucleon within a nucleus.  Check if it is
    // energetically possible to produce pions/kaons.  In not, assume nuclear excitation
    // occurs and input particle is degraded in energy. No other particles are produced.
    // If reaction is possible, find the correct number of pions/protons/neutrons
    // produced using an interpolation to multiplicity data.  Replace some pions or
    // protons/neutrons by kaons or strange baryons according to the average
    // multiplicity per Inelastic reaction.
    //
    const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass();
    const G4double etOriginal = originalIncident->GetTotalEnergy();
    const G4double pOriginal = originalIncident->GetTotalMomentum();
    const G4double targetMass = targetParticle.GetMass();
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
    G4int nt, np, nm, nz;
    const G4double c = 1.25;    
    const G4double b[] = { 0.70, 0.70 };
    if( first )       // compute normalization constants, this will only be Done once
    {
      first = false;
      G4int i;
      for( i=0; i<numMul; ++i )protmul[i] = 0.0;
      for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
      G4int counter = -1;
      for( np=0; np<(numSec/3); ++np )
      {
        for( nm=G4std::max(0,np-1); nm<=(np+1); ++nm )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nm+nz;
              if( (nt>0) && (nt<=numSec) )
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
      for( np=0; np<numSec/3; ++np )
      {
        for( nm=np; nm<=(np+2); ++nm )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nm+nz;
              if( (nt>0) && (nt<=numSec) )
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
    G4int iplab = G4std::min( 9.0, pOriginal/GeV*5.0 );
    if( (pOriginal <= 2.0*GeV) && (G4UniformRand() < cech[iplab]) )
    {
      np = nm = nz = nt = 0;
      iplab = G4std::min( 19.0, pOriginal/GeV*10.0 );
      const G4double cnk0[] = {0.17,0.18,0.17,0.24,0.26,0.20,0.22,0.21,0.34,0.45,
                               0.58,0.55,0.36,0.29,0.29,0.32,0.32,0.33,0.33,0.33};
      if( G4UniformRand() <= cnk0[iplab] )
      {
        quasiElastic = true;
        if( targetParticle.GetDefinition() == aProton )
        {
          currentParticle.SetDefinitionAndUpdateE( aKaonZL );
          incidentHasChanged = true;
          targetParticle.SetDefinitionAndUpdateE( aNeutron );
          targetHasChanged = true;
        }
      }
      else  // random number > cnk0[iplab]
      {
        G4double ran = G4UniformRand();
        if( ran < 0.25 )         // k- p --> pi- s+
        {
          if( targetParticle.GetDefinition() == aProton )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
        }
        else if( ran < 0.50 )  // k- p --> pi0 s0  or  k- n --> pi- s0
        {
          if( targetParticle.GetDefinition() == aNeutron )
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
          else
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
          targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
        else if( ran < 0.75 )  // k- p --> pi+ s-  or  k- n --> pi0 s-
        {
          if( targetParticle.GetDefinition() == aNeutron )
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
          else
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
          targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
        else                   // k- p --> pi0 L  or  k- n --> pi- L
        {
          if( targetParticle.GetDefinition() == aNeutron )
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
          else
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
          targetParticle.SetDefinitionAndUpdateE( aLambda );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
      }
    }
    else  // (pOriginal > 2.0*GeV) || (random number >= cech[iplab])
    {
      if( availableEnergy < aPiPlus->GetPDGMass() )
      {               // not energetically possible to produce pion(s)
        quasiElastic = true;
        return;
      }
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, test, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        G4int counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nm=G4std::max(0,np-1); nm<=(np+1) && ran>=excs; ++nm )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( nt > 0 )
                {
                  test = exp( G4std::min( expxu, G4std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
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
        if( np == nm )
        {
          if( G4UniformRand() >= 0.75 )
          {
            currentParticle.SetDefinitionAndUpdateE( aKaonZL );
            targetParticle.SetDefinitionAndUpdateE( aNeutron );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
        }
        else if( np == nm+1 )
        {
          targetParticle.SetDefinitionAndUpdateE( aNeutron );
          targetHasChanged = true;
        }
        else
        {
          currentParticle.SetDefinitionAndUpdateE( aKaonZL );
          incidentHasChanged = true;
        }
      }
      else   // target must be a neutron
      {
        G4int counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nm=np; nm<=(np+2) && ran>=excs; ++nm )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( (nt>=1) && (nt<=numSec) )
                {
                  test = exp( G4std::min( expxu, G4std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
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
        if( np == nm-1 )
        {
          if( G4UniformRand() < 0.5 )
          {
            targetParticle.SetDefinitionAndUpdateE( aProton );
            targetHasChanged = true;
          }
          else
          {
            currentParticle.SetDefinitionAndUpdateE( aKaonZL );
            incidentHasChanged = true;
          }
        }
        else if( np != nm )
        {
          currentParticle.SetDefinitionAndUpdateE( aKaonZL );
          incidentHasChanged = true;
        }
      }
      if( G4UniformRand() >= 0.5 )
      {
        if( (currentParticle.GetDefinition() == aKaonMinus &&
             targetParticle.GetDefinition() == aNeutron )     ||
            (currentParticle.GetDefinition() == aKaonZL &&
             targetParticle.GetDefinition() == aProton ) )
        {
          ran = G4UniformRand();
          if( ran < 0.68 )
          {
            if( targetParticle.GetDefinition() == aProton )
            {
              currentParticle.SetDefinitionAndUpdateE( aPiPlus );
              targetParticle.SetDefinitionAndUpdateE( aLambda );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
            else
            {
              currentParticle.SetDefinitionAndUpdateE( aPiMinus );
              targetParticle.SetDefinitionAndUpdateE( aLambda );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
          }
          else if( ran < 0.84 )
          {
            if( targetParticle.GetDefinition() == aProton )
            {
              currentParticle.SetDefinitionAndUpdateE( aPiZero );
              targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
            else
            {
              currentParticle.SetDefinitionAndUpdateE( aPiMinus );
              targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
          }
          else
          {
            if( targetParticle.GetDefinition() == aProton )
            {
              currentParticle.SetDefinitionAndUpdateE( aPiPlus );
              targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
            else
            {
              currentParticle.SetDefinitionAndUpdateE( aPiZero );
              targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
              incidentHasChanged = true;
              targetHasChanged = true;
            }
          }
        }
        else  // ( current != aKaonMinus || target != aNeutron ) &&
              // ( current != aKaonZL    || target != aProton  )
        {
          ran = G4UniformRand();
          if( ran < 0.67 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aLambda );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
          else if( ran < 0.78 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiMinus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaPlus );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
          else if( ran < 0.89 )
          {
            currentParticle.SetDefinitionAndUpdateE( aPiZero );
            targetParticle.SetDefinitionAndUpdateE( aSigmaZero );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
          else
          {
            currentParticle.SetDefinitionAndUpdateE( aPiPlus );
            targetParticle.SetDefinitionAndUpdateE( aSigmaMinus );
            incidentHasChanged = true;
            targetHasChanged = true;
          }
        }
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
 


