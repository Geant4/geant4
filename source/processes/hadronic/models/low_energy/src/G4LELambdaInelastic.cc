// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LELambdaInelastic.cc,v 1.1 1999-01-07 16:12:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Lambda Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 27-Mar-1997
 // Modified by J.L.Chuma 30-Apr-97: added originalTarget for CalculateMomenta
 
#include "G4LELambdaInelastic.hh"
#include "Randomize.hh"

 G4VParticleChange *
  G4LELambdaInelastic::ApplyYourself( const G4Track &aTrack,
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
      G4cout << "G4LELambdaInelastic::ApplyYourself called" << endl;
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
  G4LELambdaInelastic::Cascade(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int& vecLen,
   const G4DynamicParticle *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASL0 by H. Fesefeldt (13-Sep-1987)
    //
    // Lambda undergoes interaction with nucleon within a nucleus.  Check if it is
    // energetically possible to produce pions/kaons.  In not, assume nuclear excitation
    // occurs and input particle is degraded in energy. No other particles are produced.
    // If reaction is possible, find the correct number of pions/protons/neutrons
    // produced using an interpolation to multiplicity data.  Replace some pions or
    // protons/neutrons by kaons or strange baryons according to the average
    // multiplicity per Inelastic reaction.
    //
    const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass()/MeV;
    const G4double etOriginal = originalIncident->GetTotalEnergy()/MeV;
    const G4double targetMass = targetParticle.GetMass()/MeV;
    G4double centerofmassEnergy = sqrt( mOriginal*mOriginal +
                                        targetMass*targetMass +
                                        2.0*targetMass*etOriginal );
    G4double availableEnergy = centerofmassEnergy-(targetMass+mOriginal);
    if( availableEnergy <= G4PionPlus::PionPlus()->GetPDGMass()/MeV )
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
    G4double test;
    const G4double c = 1.25;    
    const G4double b[] = { 0.70, 0.35 };
    if( first ) {      // compute normalization constants, this will only be Done once
      first = false;
      G4int i;
      for( i=0; i<numMul; ++i )protmul[i] = 0.0;
      for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
      counter = -1;
      for( np=0; np<(numSec/3); ++np ) {
        for( nm=max(0,np-2); nm<=(np+1); ++nm ) {
          for( nz=0; nz<numSec/3; ++nz ) {
            if( ++counter < numMul ) {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec ) {
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
      for( np=0; np<numSec/3; ++np ) {
        for( nm=max(0,np-1); nm<=(np+2); ++nm ) {
          for( nz=0; nz<numSec/3; ++nz ) {
            if( ++counter < numMul ) {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec ) {
                neutmul[counter] = Pmltpc(np,nm,nz,nt,b[1],c);
                neutnorm[nt-1] += neutmul[counter];
              }
            }
          }
        }
      }
      for( i=0; i<numSec; ++i ) {
        if( protnorm[i] > 0.0 )protnorm[i] = 1.0/protnorm[i];
        if( neutnorm[i] > 0.0 )neutnorm[i] = 1.0/neutnorm[i];
      }
    }   // end of initialization
    
    const G4double expxu = 82.;           // upper bound for arg. of exp
    const G4double expxl = -expxu;        // lower bound for arg. of exp
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition *aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition *aSigmaZero = G4SigmaZero::SigmaZero();
    
    // energetically possible to produce pion(s)  -->  inelastic scattering
    //                                   otherwise quasi-elastic scattering
    
    G4double n, anpn;
    GetNormalizationConstant( availableEnergy, n, anpn );
    G4double ran = G4UniformRand();
    G4double dum, excs = 0.0;
    if( targetParticle.GetDefinition() == aProton ) {
      counter = -1;
      for( np=0; np<numSec/3 && ran>=excs; ++np ) {
        for( nm=max(0,np-2); nm<=(np+1) && ran>=excs; ++nm ) {
          for( nz=0; nz<numSec/3 && ran>=excs; ++nz ) {
            if( ++counter < numMul ) {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec ) {
                test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
                if( fabs(dum) < 1.0 ) {
                  if( test >= 1.0e-10 )excs += dum*test;
                } else {
                  excs += dum*test;
                }
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
      G4int ncht = max( 1, np-nm );
      switch( ncht ) {
       case 1:
         currentParticle.SetDefinitionAndUpdateE( aSigmaPlus );
         incidentHasChanged = true;
         break;
       case 2:
         if( G4UniformRand() < 0.5 ) {
           if( G4UniformRand() < 0.5 ) {
             currentParticle.SetDefinitionAndUpdateE( aSigmaZero );
             incidentHasChanged = true;
           }
         } else {
           currentParticle.SetDefinitionAndUpdateE( aSigmaPlus );
           incidentHasChanged = true;
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           targetHasChanged = true;
         }
         break;
       case 3:
         if( G4UniformRand() < 0.5 ) {
           if( G4UniformRand() < 0.5 ) {
             currentParticle.SetDefinitionAndUpdateE( aSigmaZero );
             incidentHasChanged = true;
           }
           targetParticle.SetDefinitionAndUpdateE( aNeutron );
           targetHasChanged = true;
         } else {
           currentParticle.SetDefinitionAndUpdateE( aSigmaMinus );
           incidentHasChanged = true;
         }
         break;
       default:
         currentParticle.SetDefinitionAndUpdateE( aSigmaMinus );
         incidentHasChanged = true;
         targetParticle.SetDefinitionAndUpdateE( aNeutron );
         targetHasChanged = true;
         break;
      }
    }
    else  // target must be a neutron
    {
      counter = -1;
      for( np=0; np<numSec/3 && ran>=excs; ++np ) {
        for( nm=max(0,np-1); nm<=(np+2) && ran>=excs; ++nm ) {
          for( nz=0; nz<numSec/3 && ran>=excs; ++nz ) {
            if( ++counter < numMul ) {
              nt = np+nm+nz;
              if( nt>0 && nt<=numSec ) {
                test = exp( min( expxu, max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                dum = (pi/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
                if( fabs(dum) < 1.0 ) {
                  if( test >= 1.0e-10 )excs += dum*test;
                } else {
                  excs += dum*test;
                }
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
      G4int ncht = max( 1, np-nm+3 );
      switch( ncht ) {
       case 1:
         currentParticle.SetDefinitionAndUpdateE( aSigmaPlus );
         incidentHasChanged = true;
         targetParticle.SetDefinitionAndUpdateE( aProton );
         targetHasChanged = true;
         break;
       case 2:
         if( G4UniformRand() < 0.5 ) {
           if( G4UniformRand() < 0.5 ) {
             currentParticle.SetDefinitionAndUpdateE( aSigmaZero );
             incidentHasChanged = true;
           }
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
         } else {
           currentParticle.SetDefinitionAndUpdateE( aSigmaPlus );
           incidentHasChanged = true;
         }
         break;
       case 3:
         if( G4UniformRand() < 0.5 ) {
           if( G4UniformRand() < 0.5 ) {
             currentParticle.SetDefinitionAndUpdateE( aSigmaZero );
             incidentHasChanged = true;
           }
         } else {
           currentParticle.SetDefinitionAndUpdateE( aSigmaMinus );
           incidentHasChanged = true;
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
         }
         break;
       default:
         currentParticle.SetDefinitionAndUpdateE( aSigmaMinus );
         incidentHasChanged = true;
         break;
      }
    }
    SetUpPions( np, nm, nz, vec, vecLen );
    return;
  }

 /* end of file */
 
