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
//
 
#include "G4RPGAntiProtonInelastic.hh"
#include "G4Exp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4HadFinalState*
G4RPGAntiProtonInelastic::ApplyYourself( const G4HadProjectile &aTrack,
                                          G4Nucleus &targetNucleus )
{ 
  const G4HadProjectile* originalIncident = &aTrack;

  if (originalIncident->GetKineticEnergy() <= 0.1*MeV) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

    //
    // create the target particle
    //
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    
    if( verboseLevel > 1 )
    {
      const G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4RPGAntiProtonInelastic::ApplyYourself called" << G4endl;
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
      G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << G4endl;
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
    modifiedOriginal.SetKineticEnergy( ek*MeV );
    et = ek + amas;
    p = std::sqrt( std::abs((et-amas)*(et+amas)) );
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
    G4FastVector<G4ReactionProduct,GHADLISTSIZE> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    const G4double cutOff = 0.1;
    const G4double anni = std::min( 1.3*originalIncident->GetTotalMomentum()/GeV, 0.4 );
    
    if( (currentParticle.GetKineticEnergy()/MeV > cutOff) ||
        (G4UniformRand() > anni) )
      Cascade( vec, vecLen,
               originalIncident, currentParticle, targetParticle,
               incidentHasChanged, targetHasChanged, quasiElastic );
    else
      quasiElastic = true;
    
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


void G4RPGAntiProtonInelastic::Cascade(
   G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
   G4int& vecLen,
   const G4HadProjectile *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
{
  // Derived from H. Fesefeldt's original FORTRAN code CASPB
  // AntiProton undergoes interaction with nucleon within a nucleus.  Check if it is
  // energetically possible to produce pions/kaons.  In not, assume nuclear excitation
  // occurs and input particle is degraded in energy. No other particles are produced.
  // If reaction is possible, find the correct number of pions/protons/neutrons
  // produced using an interpolation to multiplicity data.  Replace some pions or
  // protons/neutrons by kaons or strange baryons according to the average
  // multiplicity per Inelastic reaction.

  const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass()/MeV;
  const G4double etOriginal = originalIncident->GetTotalEnergy()/MeV;
  const G4double pOriginal = originalIncident->GetTotalMomentum()/MeV;
  const G4double targetMass = targetParticle.GetMass()/MeV;
  G4double centerofmassEnergy = std::sqrt( mOriginal*mOriginal +
                                      targetMass*targetMass +
                                      2.0*targetMass*etOriginal );
  G4double availableEnergy = centerofmassEnergy-(targetMass+mOriginal);
    
  static G4ThreadLocal G4bool first = true;
  const G4int numMul = 1200;
  const G4int numMulA = 400;
  const G4int numSec = 60;
  static G4ThreadLocal G4double protmul[numMul], protnorm[numSec]; // proton constants
  static G4ThreadLocal G4double neutmul[numMul], neutnorm[numSec]; // neutron constants
  static G4ThreadLocal G4double protmulA[numMulA], protnormA[numSec]; // proton constants
  static G4ThreadLocal G4double neutmulA[numMulA], neutnormA[numSec]; // neutron constants
  // np = number of pi+, nneg = number of pi-, nz = number of pi0
  G4int counter, nt=0, np=0, nneg=0, nz=0;
  G4double test;
  const G4double c = 1.25;    
  const G4double b[] = { 0.7, 0.7 };

  if (first) {       // compute normalization constants, this will only be Done once
    first = false;
    G4int i;
    for( i=0; i<numMul; ++i )protmul[i] = 0.0;
    for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
    counter = -1;
    for (np=0; np<(numSec/3); ++np ) {
      for (nneg=std::max(0,np-1); nneg<=(np+1); ++nneg ) {
        for (nz=0; nz<numSec/3; ++nz ) {
          if (++counter < numMul ) {
            nt = np+nneg+nz;
            if (nt>0 && nt<=numSec ) {
              protmul[counter] = Pmltpc(np,nneg,nz,nt,b[0],c);
              protnorm[nt-1] += protmul[counter];
            }
          }
        }
      }
    }

    for (i=0; i<numMul; ++i )neutmul[i] = 0.0;
    for (i=0; i<numSec; ++i )neutnorm[i] = 0.0;
    counter = -1;
    for (np=0; np<numSec/3; ++np ) {
        for( nneg=np; nneg<=(np+2); ++nneg )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nneg+nz;
              if( (nt>0) && (nt<=numSec) )
              {
                neutmul[counter] = Pmltpc(np,nneg,nz,nt,b[1],c);
                neutnorm[nt-1] += neutmul[counter];
              }
            }
          }
        }
    }

    for (i=0; i<numSec; ++i ) {
      if( protnorm[i] > 0.0 )protnorm[i] = 1.0/protnorm[i];
      if( neutnorm[i] > 0.0 )neutnorm[i] = 1.0/neutnorm[i];
    }

    // do the same for annihilation channels

    for (i=0; i<numMulA; ++i )protmulA[i] = 0.0;
    for (i=0; i<numSec; ++i )protnormA[i] = 0.0;
    counter = -1;

    for (np=0; np<(numSec/3); ++np ) {
      nneg = np;
      for (nz=0; nz<numSec/3; ++nz ) {
        if ( ++counter < numMulA ) {
          nt = np+nneg+nz;
          if ( nt>1 && nt<=numSec ) {
            protmulA[counter] = Pmltpc(np,nneg,nz,nt,b[0],c);
            protnormA[nt-1] += protmulA[counter];
          }
        }
      }
    }

    for( i=0; i<numMulA; ++i )neutmulA[i] = 0.0;
    for( i=0; i<numSec; ++i )neutnormA[i] = 0.0;
    counter = -1;

    for (np=0; np<numSec/3; ++np ) {
      nneg = np+1;
      for (nz=0; nz<numSec/3; ++nz ) {
        if (++counter < numMulA ) {
          nt = np+nneg+nz;
          if ( (nt>1) && (nt<=numSec) ) {
            neutmulA[counter] = Pmltpc(np,nneg,nz,nt,b[1],c);
            neutnormA[nt-1] += neutmulA[counter];
          }
        }
      }
    }
    for (i=0; i<numSec; ++i) {
      if( protnormA[i] > 0.0 )protnormA[i] = 1.0/protnormA[i];
      if( neutnormA[i] > 0.0 )neutnormA[i] = 1.0/neutnormA[i];
    }
  }   // end of initialization

  // initialization is OK

  const G4double expxu = 82.;           // upper bound for arg. of exp
  const G4double expxl = -expxu;        // lower bound for arg. of exp
    
  G4ParticleDefinition *aNeutron      = G4Neutron::Neutron();
  G4ParticleDefinition *aProton       = G4Proton::Proton();
  G4ParticleDefinition *anAntiNeutron = G4AntiNeutron::AntiNeutron();
  G4ParticleDefinition *aPiPlus       = G4PionPlus::PionPlus();
    
  const G4double anhl[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.90,
                           0.6,0.52,0.47,0.44,0.41,0.39,0.37,0.35,0.34,0.24,
                           0.19,0.15,0.12,0.10,0.09,0.07,0.06,0.05,0.0};
    
  G4int iplab = G4int( pOriginal/GeV*10.0 );
  if( iplab >  9 )iplab = G4int( pOriginal/GeV ) + 9;
  if( iplab > 18 )iplab = G4int( pOriginal/GeV/10.0 ) + 18;
  if( iplab > 27 )iplab = 28;
  if (G4UniformRand() > anhl[iplab] ) {
    if (availableEnergy <= aPiPlus->GetPDGMass()/MeV ) {
      quasiElastic = true;
      return;
    }
    G4int ieab = static_cast<G4int>(availableEnergy*5.0/GeV);
    const G4double supp[] = {0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98};
    G4double w0, wp, wt, wm;
    if ( (availableEnergy < 2.0*GeV) && (G4UniformRand() >= supp[ieab]) ) {
      //
      // suppress high multiplicity events at low momentum
      // only one pion will be produced
      //
      np = nneg = nz = 0;
      if ( targetParticle.GetDefinition() == aProton ) {
        test = G4Exp( std::min( expxu, std::max( expxl, -(1.0+b[1])*(1.0+b[1])/(2.0*c*c) ) ) );
        w0 = test;
        wp = test;        
        test = G4Exp( std::min( expxu, std::max( expxl, -(-1.0+b[1])*(-1.0+b[1])/(2.0*c*c) ) ) );
        wm = test;
        wt = w0+wp+wm;
        wp += w0;
        G4double ran = G4UniformRand();
        if( ran < w0/wt )
          nz = 1;
        else if( ran < wp/wt )
          np = 1;
        else
          nneg = 1;
      } else { // target is a neutron
        test = G4Exp( std::min( expxu, std::max( expxl, -(1.0+b[0])*(1.0+b[0])/(2.0*c*c) ) ) );
        w0 = test;
        test = G4Exp( std::min( expxu, std::max( expxl, -(-1.0+b[0])*(-1.0+b[0])/(2.0*c*c) ) ) );
        wm = test;
        G4double ran = G4UniformRand();
        if( ran < w0/(w0+wm) )
          nz = 1;
        else
          nneg = 1;
      }
    } else {  // (availableEnergy >= 2.0*GeV) || (random number < supp[ieab])
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if (targetParticle.GetDefinition() == aProton ) {
        counter = -1;
        for (np=0; np<numSec/3 && ran>=excs; ++np ) {
          for (nneg=std::max(0,np-1); nneg<=(np+1) && ran>=excs; ++nneg ) {
            for (nz=0; nz<numSec/3 && ran>=excs; ++nz ) {
              if (++counter < numMul ) {
                nt = np+nneg+nz;
                if ((nt>0) && (nt<=numSec) ) {
                  test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
                  if (std::fabs(dum) < 1.0 ) {
                    if( test >= 1.0e-10 )excs += dum*test;
                  }
                    else
                      excs += dum*test;
                }
              }
            }
          }
        }
        if (ran >= excs )  { // 3 previous loops continued to the end
          quasiElastic = true;
          return;
        }
      } else { // target must be a neutron
        counter = -1;
        for ( np=0; np<numSec/3 && ran>=excs; ++np ) {
          for ( nneg=np; nneg<=(np+2) && ran>=excs; ++nneg ) {
            for ( nz=0; nz<numSec/3 && ran>=excs; ++nz ) {
              if ( ++counter < numMul ) {
                nt = np+nneg+nz;
                if ( (nt>0) && (nt<=numSec) ) {
                  test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
                  if ( std::fabs(dum) < 1.0 )
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
        if( ran >= excs ) { // 3 previous loops continued to the end
          quasiElastic = true;
          return;
        }
      }
      np--; nneg--; nz--;

    } // end if availableEnergy

    if (targetParticle.GetDefinition() == aProton ) {
      switch( np-nneg )
      {
       case 0:
         if( G4UniformRand() < 0.33 )
         {
           currentParticle.SetDefinitionAndUpdateE( anAntiNeutron );
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
         currentParticle.SetDefinitionAndUpdateE( anAntiNeutron );
         incidentHasChanged = true;
         break;
      }

    } else { // target must be a neutron
      switch( np-nneg )
      {
       case -1:
           if( G4UniformRand() < 0.5 )
           {
             targetParticle.SetDefinitionAndUpdateE( aProton );
             targetHasChanged = true;
           }
           else
           {
             currentParticle.SetDefinitionAndUpdateE( anAntiNeutron );
             incidentHasChanged = true;
           }
           break;
       case 0:
           break;
       default:
           currentParticle.SetDefinitionAndUpdateE( anAntiNeutron );
           targetParticle.SetDefinitionAndUpdateE( aProton );
           incidentHasChanged = true;
           targetHasChanged = true;
           break;
      }
    }

  } else {  //  random number <= anhl[iplab]
    if (centerofmassEnergy <= 2*aPiPlus->GetPDGMass()/MeV ) {
      quasiElastic = true;
      return;
    }

    // annihilation channels

    G4double n, anpn;
    GetNormalizationConstant( -centerofmassEnergy, n, anpn );
    G4double ran = G4UniformRand();
    G4double dum, excs = 0.0;
    if (targetParticle.GetDefinition() == aProton ) {
        counter = -1;
        for( np=0; (np<numSec/3) && (ran>=excs); ++np )
        {
          nneg = np;
          for( nz=0; (nz<numSec/3) && (ran>=excs); ++nz )
          {
            if( ++counter < numMulA )
            {
              nt = np+nneg+nz;
              if( (nt>1) && (nt<=numSec) )
              {
                test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*protmulA[counter]*protnormA[nt-1]/(2.0*n*n);
                if( std::abs(dum) < 1.0 )
                {
                  if( test >= 1.0e-10 )excs += dum*test;
                }
                else
                  excs += dum*test;
              }
            }
          }
        }
    } else { // target must be a neutron
      counter = -1;
        for( np=0; (np<numSec/3) && (ran>=excs); ++np )
        {
          nneg = np+1;
          for( nz=0; (nz<numSec/3) && (ran>=excs); ++nz )
          {
            if( ++counter < numMulA )
            {
              nt = np+nneg+nz;
              if( (nt>1) && (nt<=numSec) )
              {
                test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*neutmulA[counter]*neutnormA[nt-1]/(2.0*n*n);
                if( std::fabs(dum) < 1.0 )
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
    if ( ran >= excs ) { // 3 previous loops continued to the end
      quasiElastic = true;
      return;
    }
    np--; nz--;
    currentParticle.SetMass( 0.0 );
    targetParticle.SetMass( 0.0 );
  }

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl;
  while(np + nneg + nz < 3) {   /* Loop checking, 01.09.2015, D.Wright */
    nz++;
    loop++;
    if (loop > 1000) {
      G4Exception("G4RPGAntiProtonInelastic::Cascade()", "HAD_RPG_100", JustWarning, ed);
      break;
    }
  }

  SetUpPions( np, nneg, nz, vec, vecLen );
  return;
}

 /* end of file */
 
