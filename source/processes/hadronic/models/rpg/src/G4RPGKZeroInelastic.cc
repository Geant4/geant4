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
// $Id: G4RPGKZeroInelastic.cc 94214 2015-11-09 08:18:05Z gcosmo $
//
 
#include "G4RPGKZeroInelastic.hh"
#include "G4Exp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
 
G4HadFinalState*
G4RPGKZeroInelastic::ApplyYourself( const G4HadProjectile &aTrack,
                                        G4Nucleus &targetNucleus )
{
  const G4HadProjectile *originalIncident = &aTrack;

  // create the target particle

  G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    
  if( verboseLevel > 1 )
  {
    const G4Material *targetMaterial = aTrack.GetMaterial();
    G4cout << "G4RPGKZeroInelastic::ApplyYourself called" << G4endl;
    G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
    G4cout << "target material = " << targetMaterial->GetName() << ", ";
    G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << G4endl;
  }

  // Fermi motion and evaporation
  // As of Geant3, the Fermi energy calculation had not been Done

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
 
void G4RPGKZeroInelastic::Cascade(
   G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
   G4int& vecLen,
   const G4HadProjectile *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
{
  // Derived from H. Fesefeldt's original FORTRAN code CASK0
  //
  // K0Short undergoes interaction with nucleon within a nucleus.  Check if it is
  // energetically possible to produce pions/kaons.  In not, assume nuclear excitation
  // occurs and input particle is degraded in energy. No other particles are produced.
  // If reaction is possible, find the correct number of pions/protons/neutrons
  // produced using an interpolation to multiplicity data.  Replace some pions or
  // protons/neutrons by kaons or strange baryons according to the average
  // multiplicity per Inelastic reaction.

  const G4double mOriginal = originalIncident->GetDefinition()->GetPDGMass()/MeV;
  const G4double etOriginal = originalIncident->GetTotalEnergy()/MeV;
  const G4double targetMass = targetParticle.GetMass()/MeV;
  G4double centerofmassEnergy = std::sqrt( mOriginal*mOriginal +
                                      targetMass*targetMass +
                                      2.0*targetMass*etOriginal );
  G4double availableEnergy = centerofmassEnergy-(targetMass+mOriginal);
  if( availableEnergy <= G4PionPlus::PionPlus()->GetPDGMass()/MeV )
  {
    quasiElastic = true;
    return;
  }
    static G4ThreadLocal G4bool first = true;
    const G4int numMul = 1200;
    const G4int numSec = 60;
    static G4ThreadLocal G4double protmul[numMul], protnorm[numSec]; // proton constants
    static G4ThreadLocal G4double neutmul[numMul], neutnorm[numSec]; // neutron constants
    // np = number of pi+, nneg = number of pi-, nz = number of pi0
    G4int counter, nt=0, np=0, nneg=0, nz=0;
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
        for( nneg=std::max(0,np-1); nneg<=(np+1); ++nneg )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nneg+nz;
              if( nt>0 && nt<=numSec )
              {
                protmul[counter] = Pmltpc(np,nneg,nz,nt,b[0],c);
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
        for( nneg=std::max(0,np-2); nneg<=np; ++nneg )
        {
          for( nz=0; nz<numSec/3; ++nz )
          {
            if( ++counter < numMul )
            {
              nt = np+nneg+nz;
              if( nt>0 && nt<=numSec )
              {
                neutmul[counter] = Pmltpc(np,nneg,nz,nt,b[1],c);
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
    G4ParticleDefinition *aKaonZL = G4KaonZeroLong::KaonZeroLong();
    G4ParticleDefinition *aKaonZS = G4KaonZeroShort::KaonZeroShort();
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4int ieab = static_cast<G4int>(5.0*availableEnergy*MeV/GeV);
    const G4double supp[] = {0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98};
    G4double test, w0, wp, wt, wm;
    if( (availableEnergy*MeV/GeV < 2.0) && (G4UniformRand() >= supp[ieab]) )
    {
      //
      // suppress high multiplicity events at low momentum
      // only one pion will be produced
      //
      nneg = np = nz = 0;
      if( targetParticle.GetDefinition() == aNeutron )
      {
        test = G4Exp( std::min( expxu, std::max( expxl, -(1.0+b[0])*(1.0+b[0])/(2.0*c*c) ) ) );
        w0 = test/2.0;
        test = G4Exp( std::min( expxu, std::max( expxl, -(-1.0+b[0])*(1.0+b[0])/(2.0*c*c) ) ) );
        wm = test*1.5;
        if( G4UniformRand() < w0/(w0+wm) )
          nz = 1;
        else
          nneg = 1;
      }
      else  // target is a proton
      {
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
      }
    }
    else //  (availableEnergy*MeV/GeV >= 2.0) || (G4UniformRand() < supp[ieab])
    {
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nneg=std::max(0,np-1); nneg<=(np+1) && ran>=excs; ++nneg )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nneg+nz;
                if( nt>0 && nt<=numSec )
                {
                  test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
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
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nneg--; nz--;
      }
      else  // target must be a neutron
      {
        counter = -1;
        for( np=0; np<numSec/3 && ran>=excs; ++np )
        {
          for( nneg=std::max(0,np-2); nneg<=np && ran>=excs; ++nneg )
          {
            for( nz=0; nz<numSec/3 && ran>=excs; ++nz )
            {
              if( ++counter < numMul )
              {
                nt = np+nneg+nz;
                if( nt>0 && nt<=numSec )
                {
                  test = G4Exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
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
        if( ran >= excs )  // 3 previous loops continued to the end
        {
          quasiElastic = true;
          return;
        }
        np--; nneg--; nz--;
      }
    }
    if( targetParticle.GetDefinition() == aProton )
    {
      switch( np-nneg )
      {
       case 0:
         if( G4UniformRand() < 0.25 )
         {
           currentParticle.SetDefinitionAndUpdateE( aKaonPlus );
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
         targetParticle.SetDefinitionAndUpdateE( aNeutron );
         targetHasChanged = true;
         break;
      }
    }
    else   // targetParticle is a neutron
    {
      switch( np-nneg )          // seems wrong, charge not conserved
      {
       case 1:
         if( G4UniformRand() < 0.5 )
         {
           currentParticle.SetDefinitionAndUpdateE( aKaonPlus );
           incidentHasChanged = true;
         }
         else
         {
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
         }
         break;
       case 2:
         currentParticle.SetDefinitionAndUpdateE( aKaonPlus );
         incidentHasChanged = true;
         targetParticle.SetDefinitionAndUpdateE( aProton );
         targetHasChanged = true;
         break;
       default:
         break;
      }
    }

  if (currentParticle.GetDefinition() == aKaonZS) {
    if (G4UniformRand() >= 0.5) {
      currentParticle.SetDefinitionAndUpdateE(aKaonZL);
      incidentHasChanged = true;
    }
  }

  if (targetParticle.GetDefinition() == aKaonZS) {
    if (G4UniformRand() >= 0.5) {
      targetParticle.SetDefinitionAndUpdateE(aKaonZL);
      targetHasChanged = true;
    }
  }

  SetUpPions( np, nneg, nz, vec, vecLen );
  return;
}

 /* end of file */
 
