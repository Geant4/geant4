// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEPionMinusInelastic.cc,v 1.1 1999-01-07 16:12:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: PionMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Nov-1996
 // Last modified: 27-Mar-1997
 // Modified by J.L.Chuma 30-Apr-97: added originalTarget for CalculateMomenta
 
#include "G4LEPionMinusInelastic.hh"
#include "Randomize.hh"
 
 G4VParticleChange *
  G4LEPionMinusInelastic::ApplyYourself( const G4Track &aTrack,
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
      G4cout << "G4PionMinusInelastic::ApplyYourself called" << endl;
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy() << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
      G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << endl;
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
  G4LEPionMinusInelastic::Cascade(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int& vecLen,
   const G4DynamicParticle *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASPIM by H. Fesefeldt (13-Sep-1987)
    //
    // pi-  undergoes interaction with nucleon within nucleus.
    // Check if energetically possible to produce pions/kaons.
    // If not assume nuclear excitation occurs and input particle
    // is degraded in energy.  No other particles produced.
    // If reaction is possible find correct number of pions/protons/neutrons
    // produced using an interpolation to multiplicity data.
    // Replace some pions or protons/neutrons by kaons or strange baryons
    // according to average multiplicity per inelastic reactions.
    //
    const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass();
    const G4double etOriginal = originalIncident->GetTotalEnergy();
    const G4double pOriginal = originalIncident->GetTotalMomentum();
    const G4double targetMass = targetParticle.GetMass();
    G4double centerofmassEnergy = sqrt( mOriginal*mOriginal +
                                        targetMass*targetMass +
                                        2.0*targetMass*etOriginal );
    G4double availableEnergy = centerofmassEnergy-(targetMass+mOriginal);
    if( availableEnergy <= G4PionMinus::PionMinus()->GetPDGMass() )
    {
      quasiElastic = true;
      return;
    }
    static G4bool first = true;
    const G4int numMul = 1200;
    const G4int numSec = 60;
    static G4double protmul[numMul], protnorm[numSec]; // proton constants
    static G4double neutmul[numMul], neutnorm[numSec]; // neutron constants
    // np = number of pi+, nm = number of pi-, nz = number of pi0
    G4int counter, nt=0, np=0, nm=0, nz=0;
    const G4double c = 1.25;
    const G4double b[] = { 0.70, 0.70 };
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
              if( nt > 0 )
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
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aPiZero = G4PionZero::PionZero();
    G4int ieab = availableEnergy*5.0/GeV;
    const G4double supp[] = {0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98};
    G4double test, w0, wp, wt, wm;
    if( (availableEnergy<2.0*GeV) && (G4UniformRand()>=supp[ieab]) )
    {
      // suppress high multiplicity events at low momentum
      // only one pion will be produced
      
      // charge exchange reaction is included in inelastic cross section
        
      const G4double cech[] = {1.,0.95,0.79,0.32,0.19,0.16,0.14,0.12,0.10,0.08};
      G4int iplab = min( 9.0, pOriginal/GeV*5.0 );
      if( G4UniformRand() <= cech[iplab] )
      {
        if( targetParticle.GetDefinition() == aProton )
        {
          currentParticle.SetDefinitionAndUpdateE( aPiZero );  // charge exchange
          targetParticle.SetDefinitionAndUpdateE( aNeutron );
          incidentHasChanged = true;
          targetHasChanged = true;
        }
        else quasiElastic = true;
        return;
      }
      nm = np = nz = 0;
      if( targetParticle.GetDefinition() == aProton )
      {
        test = exp( min( expxu, max( expxl, -(1.0+b[0])*(1.0+b[0])/(2.0*c*c) ) ) );
        w0 = test;
        wp = 10.0*test;        
        test = exp( min( expxu, max( expxl, -(-1.0+b[0])*(-1.0+b[0])/(2.0*c*c) ) ) );
        wm = test;
        wt = w0+wp+wm;
        wp += w0;
        G4double ran = G4UniformRand();
        if( ran < w0/wt )
          nz = 1;
        else if( ran < wp/wt )
          np = 1;
        else
          nm = 1;
      }
      else  // target is a neutron
      {
        test = exp( min( expxu, max( expxl, -(1.0+b[1])*(1.0+b[1])/(2.0*c*c) ) ) );
        w0 = test;
        test = exp( min( expxu, max( expxl, -(-1.0+b[1])*(-1.0+b[1])/(2.0*c*c) ) ) );
        wm = test;
        G4double ran = G4UniformRand();
        if( ran < w0/(w0+wm) )
          nz = 1;
        else
          nm = 1;
      }
    }
    else
    {
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        counter = -1;
        for( np=0; (np<numSec/3) && (ran>=excs); ++np )
        {
          for( nm=max(0,np-1); (nm<=(np+1)) && (ran>=excs); ++nm )
          {
            for( nz=0; (nz<numSec/3) && (ran>=excs); ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( nt > 0 )
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
      }
      else  // target must be a neutron
      {
        counter = -1;
        for( np=0; (np<numSec/3) && (ran>=excs); ++np )
        {
          for( nm=np; (nm<=(np+2)) && (ran>=excs); ++nm )
          {
            for( nz=0; (nz<numSec/3) && (ran>=excs); ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nm+nz;
                if( (nt>=1) && (nt<=numSec) )
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
      }
    }
    if( targetParticle.GetDefinition() == aProton )
    {
      switch( np-nm )
      {
       case 0:
         if( G4UniformRand() >= 0.75 )
         {
           currentParticle.SetDefinitionAndUpdateE( aPiZero );
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           incidentHasChanged = true;
           targetHasChanged = true;
         }
         break;
       case 1:
         targetParticle.SetDefinitionAndUpdateE( aNeutron );
         targetHasChanged = true;
         break;
       default:
         currentParticle.SetDefinitionAndUpdateE( aPiZero );
         incidentHasChanged = true;
         break;
      }
    }
    else
    {
      switch( np-nm )
      {
       case -1:
         if( G4UniformRand() < 0.5 )
         {
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
         } else {
           currentParticle.SetDefinitionAndUpdateE( aPiZero );
           incidentHasChanged = true;
         }
         break;
       case 0:
         break;
       default:
         currentParticle.SetDefinitionAndUpdateE( aPiZero );
         incidentHasChanged = true;
         break;
      }
    }
    SetUpPions( np, nm, nz, vec, vecLen );
    return;
  }

 /* end of file */
 
