// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEAntiSigmaPlusInelastic.cc,v 1.1 1999-01-07 16:12:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiSigmaPlus Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Feb-1997
 // Last modified: 27-Mar-1997
 // Modified by J.L.Chuma 30-Apr-97: added originalTarget for CalculateMomenta
 
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "Randomize.hh"

 G4VParticleChange *
  G4LEAntiSigmaPlusInelastic::ApplyYourself( const G4Track &aTrack,
                                             G4Nucleus &targetNucleus )
  {
    theParticleChange.Initialize( aTrack );
    
    const G4DynamicParticle *originalIncident = aTrack.GetDynamicParticle();
    if (originalIncident->GetKineticEnergy()<= 0.1*MeV) return &theParticleChange;
    //
    // create the target particle
    //
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    
    if( verboseLevel > 1 )
    {
      G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEAntiSigmaPlusInelastic::ApplyYourself called" << endl;
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
    const G4double anni = min( 1.3*currentParticle.GetTotalMomentum()/GeV, 0.4 );
    if( (currentParticle.GetKineticEnergy()/MeV > cutOff) || (G4UniformRand() > anni) )
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
  G4LEAntiSigmaPlusInelastic::Cascade(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int& vecLen,
   const G4DynamicParticle *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASASP by H. Fesefeldt (13-Sep-1987)
    //
    // AntiSigmaPlus undergoes interaction with nucleon within a nucleus.  Check if it is
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
    const G4int numMulA = 400;
    const G4int numSec = 60;
    static G4double protmul[numMul], protnorm[numSec]; // proton constants
    static G4double neutmul[numMul], neutnorm[numSec]; // neutron constants
    static G4double protmulA[numMulA], protnormA[numSec]; // proton constants
    static G4double neutmulA[numMulA], neutnormA[numSec]; // neutron constants
    // np = number of pi+, nm = number of pi-, nz = number of pi0
    G4int counter, nt=0, np=0, nm=0, nz=0;
    G4double test;
    const G4double c = 1.25;    
    const G4double b[] = { 0.7, 0.7 };
    if( first )       // compute normalization constants, this will only be Done once
    {
      first = false;
      G4int i;
      for( i=0; i<numMul; ++i )protmul[i] = 0.0;
      for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
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
      //
      // do the same for annihilation channels
      //
      for( i=0; i<numMulA; ++i )protmulA[i] = 0.0;
      for( i=0; i<numSec; ++i )protnormA[i] = 0.0;
      counter = -1;
      for( np=1; np<(numSec/3); ++np )
      {
        nm = np;
        for( nz=0; nz<numSec/3; ++nz )
        {
          if( ++counter < numMulA )
          {
            nt = np+nm+nz;
            if( nt>1 && nt<=numSec )
            {
              protmulA[counter] = Pmltpc(np,nm,nz,nt,b[0],c);
              protnormA[nt-1] += protmulA[counter];
            }
          }
        }
      }
      for( i=0; i<numMulA; ++i )neutmulA[i] = 0.0;
      for( i=0; i<numSec; ++i )neutnormA[i] = 0.0;
      counter = -1;
       for( np=0; np<numSec/3; ++np )
       {
        nm = np+1;
        for( nz=0; nz<numSec/3; ++nz )
        {
          if( ++counter < numMulA )
          {
            nt = np+nm+nz;
            if( nt>1 && nt<=numSec )
            {
              neutmulA[counter] = Pmltpc(np,nm,nz,nt,b[1],c);
              neutnormA[nt-1] += neutmulA[counter];
            }
          }
        }
      }
      for( i=0; i<numSec; ++i )
      {
        if( protnormA[i] > 0.0 )protnormA[i] = 1.0/protnormA[i];
        if( neutnormA[i] > 0.0 )neutnormA[i] = 1.0/neutnormA[i];
      }
    }   // end of initialization
    
    const G4double expxu = 82.;           // upper bound for arg. of exp
    const G4double expxl = -expxu;        // lower bound for arg. of exp
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aPiPlus = G4PionPlus::PionPlus();
    G4ParticleDefinition *anAntiLambda = G4AntiLambda::AntiLambda();
    G4ParticleDefinition *aKaonMinus = G4KaonMinus::KaonMinus();
    G4ParticleDefinition *aKaonPlus = G4KaonPlus::KaonPlus();
    G4ParticleDefinition *aKaonZL = G4KaonZeroLong::KaonZeroLong();
    G4ParticleDefinition *anAntiSigmaZero = G4AntiSigmaZero::AntiSigmaZero();
    const G4double anhl[] = {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.97,0.88,
                             0.85,0.81,0.75,0.64,0.64,0.55,0.55,0.45,0.47,0.40,
                             0.39,0.36,0.33,0.10,0.01};
    G4int iplab = G4int( pOriginal/GeV*10.0 );
    if( iplab >  9 )iplab = G4int( (pOriginal/GeV- 1.0)*5.0  ) + 10;
    if( iplab > 14 )iplab = G4int(  pOriginal/GeV- 2.0       ) + 15;
    if( iplab > 22 )iplab = G4int( (pOriginal/GeV-10.0)/10.0 ) + 23;
    if( iplab > 24 )iplab = 24;
    if( G4UniformRand() > anhl[iplab] )
    {
      if( availableEnergy <= aPiPlus->GetPDGMass()/MeV )
      {
        quasiElastic = true;
        return;
      }
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
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
                if( (nt>0) && (nt<=numSec) )
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
        G4int ncht = min( 3, max( 1, np-nm+2 ) );
        switch( ncht )
        {
         case 1:
           if( G4UniformRand() < 0.5 )
             currentParticle.SetDefinitionAndUpdateE( anAntiLambda );
           else
             currentParticle.SetDefinitionAndUpdateE( anAntiSigmaZero );
           incidentHasChanged = true;
           break;
         case 2:
           if( G4UniformRand() >= 0.5 )
           {
             if( G4UniformRand() < 0.5 )
               currentParticle.SetDefinitionAndUpdateE( anAntiLambda );
             else
               currentParticle.SetDefinitionAndUpdateE( anAntiSigmaZero );
             incidentHasChanged = true;
           }             
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           targetHasChanged = true;
           break;
         case 3:
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           targetHasChanged = true;
           break;
        }
      }
      else  // target must be a neutron
      {
        counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nm=np; nm<=(np+2) && ran>=excs; ++nm )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( (nt>0) && (nt<=numSec) )
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
        G4int ncht = min( 3, max( 1, np-nm+3 ) );
        switch( ncht )
        {
         case 1:
           if( G4UniformRand() < 0.5 )
             currentParticle.SetDefinitionAndUpdateE( anAntiLambda );
           else
             currentParticle.SetDefinitionAndUpdateE( anAntiSigmaZero );
           incidentHasChanged = true;
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
           break;
         case 2:
           if( G4UniformRand() < 0.5 )
           {
             if( G4UniformRand() < 0.5 )
             {
               currentParticle.SetDefinitionAndUpdateE( anAntiLambda );
               incidentHasChanged = true;
             }
             else
             {
               targetParticle.SetDefinitionAndUpdateE( aProton );
               targetHasChanged = true;
             }
           }
           else
           {
             if( G4UniformRand() < 0.5 )
             {
               currentParticle.SetDefinitionAndUpdateE( anAntiSigmaZero );
               incidentHasChanged = true;
             }
             else
             {
               targetParticle.SetDefinitionAndUpdateE( aProton );
               targetHasChanged = true;
             }
           }
           break;
         case 3:
           break;
        }
      }
    }
    else  // random number <= anhl[iplab]
    {
      if( centerofmassEnergy <= aPiPlus->GetPDGMass()/MeV+aKaonPlus->GetPDGMass()/MeV )
      {
        quasiElastic = true;
        return;
      }
      G4double n, anpn;
      GetNormalizationConstant( -centerofmassEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        counter = -1;
        for( np=1; np<numSec/3 && ran>=excs; ++np )
        {
          nm = np;
          for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
          {
            if( ++counter < numMulA )
            {
              nt = np+nm+nz;
              if( nt>1 && nt<=numSec )
              {
                test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*protmulA[counter]*protnormA[nt-1]/(2.0*n*n);
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
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nz--;
      }
      else  // target must be a neutron
      {
        counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          nm = np+1;
          for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
          {
            if( ++counter < numMulA )
            {
              nt = np+nm+nz;
              if( nt>1 && nt<=numSec )
              {
                test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*neutmulA[counter]*neutnormA[nt-1]/(2.0*n*n);
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
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nz--;
      }
      if( nz > 0 )
      {
        if( nm > 0 )
        {
          if( G4UniformRand() < 0.5 )
          {
            vec.Initialize( 1 );
            G4ReactionProduct *p= new G4ReactionProduct;
            p->SetDefinition( aKaonMinus );
            (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
            vec.SetElement( vecLen++, p );
            --nm;
          }
          else
          {
            vec.Initialize( 1 );
            G4ReactionProduct *p= new G4ReactionProduct ;
            p->SetDefinition( aKaonZL );
            (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
            vec.SetElement( vecLen++, p );
            --nz;
          }
        }
        else   // nm == 0
        {
          vec.Initialize( 1 );
          G4ReactionProduct *p = new G4ReactionProduct;
          p->SetDefinition( aKaonZL );
          (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
          vec.SetElement( vecLen++, p );
          --nz;
        }
      }
      else    //  nz == 0
      {
        if( nm > 0 )
        {
          vec.Initialize( 1 );
          G4ReactionProduct *p = new G4ReactionProduct;
          p->SetDefinition( aKaonMinus );
          (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
          vec.SetElement( vecLen++, p );
          --nm;
        }
      }
      currentParticle.SetMass( 0.0 );
      targetParticle.SetMass( 0.0 );
    }
    SetUpPions( np, nm, nz, vec, vecLen );
    return;
  }

 /* end of file */
 
