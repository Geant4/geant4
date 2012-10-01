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
// Hadronic Process: Low Energy Neutron Inelastic Process
// J.L. Chuma, TRIUMF, 04-Feb-1997
 
#include <iostream>

#include "G4LENeutronInelastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

void G4LENeutronInelastic::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4LENeutronInelastic is one of the Low Energy Parameterized\n"
          << "(LEP) models used to implement inelastic neutron scattering\n"
          << "from nuclei.  It is a re-engineered version of the GHEISHA\n"
          << "code of H. Fesefeldt.  It divides the initial collision\n"
          << "products into backward- and forward-going clusters which are\n"
          << "then decayed into final state hadrons.  The model does not\n"
          << "conserve energy on an event-by-event basis.  It may be\n"
          << "applied to neutrons with initial energies between 0 and 25\n"
          << "GeV.\n";
}


G4HadFinalState*
G4LENeutronInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                    G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile *originalIncident = &aTrack;

  // Create the target particle
  G4DynamicParticle* originalTarget = targetNucleus.ReturnTargetParticle();
    
  if (verboseLevel > 1) {
    const G4Material* targetMaterial = aTrack.GetMaterial();
    G4cout << "G4LENeutronInelastic::ApplyYourself called" << G4endl;
    G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
    G4cout << "target material = " << targetMaterial->GetName() << ", ";
    G4cout << "target particle = " << originalTarget->GetDefinition()->GetParticleName()
           << G4endl;
  }
 
  G4ReactionProduct modifiedOriginal;
  modifiedOriginal = *originalIncident;
  G4ReactionProduct targetParticle;
  targetParticle = *originalTarget;
  if (originalIncident->GetKineticEnergy()/GeV < 0.01 + 2.*G4UniformRand()/9.) {
    SlowNeutron(originalIncident, modifiedOriginal, targetParticle, targetNucleus);
    if (isotopeProduction) DoIsotopeCounting(originalIncident, targetNucleus);
    delete originalTarget;
    return &theParticleChange;
  }

  // Fermi motion and evaporation
  // As of Geant3, the Fermi energy calculation had not been done
  G4double ek = originalIncident->GetKineticEnergy()/MeV;
  G4double amas = originalIncident->GetDefinition()->GetPDGMass()/MeV;
    
  G4double tkin = targetNucleus.Cinema(ek);
  ek += tkin;
  modifiedOriginal.SetKineticEnergy( ek*MeV );
  G4double et = ek + amas;
  G4double p = std::sqrt( std::abs((et-amas)*(et+amas)) );
  G4double pp = modifiedOriginal.GetMomentum().mag()/MeV;
  if (pp > 0.0) {
    G4ThreeVector momentum = modifiedOriginal.GetMomentum();
    modifiedOriginal.SetMomentum(momentum * (p/pp) );
  }

  // calculate black track energies
  tkin = targetNucleus.EvaporationEffects( ek );
  ek -= tkin;
  modifiedOriginal.SetKineticEnergy( ek*MeV );
  et = ek + amas;
  p = std::sqrt( std::abs((et-amas)*(et+amas)) );
  pp = modifiedOriginal.GetMomentum().mag()/MeV;
  if (pp > 0.0) {
    G4ThreeVector momentum = modifiedOriginal.GetMomentum();
    modifiedOriginal.SetMomentum(momentum * (p/pp) );
  }
  const G4double cutOff = 0.1;
  if (modifiedOriginal.GetKineticEnergy()/MeV <= cutOff) {
    SlowNeutron(originalIncident, modifiedOriginal, targetParticle, targetNucleus);
    if (isotopeProduction) DoIsotopeCounting(originalIncident, targetNucleus);
    delete originalTarget;
    return &theParticleChange;
  }

  G4ReactionProduct currentParticle = modifiedOriginal;
  currentParticle.SetSide(1);  // incident always goes in forward hemisphere
  targetParticle.SetSide(-1);  // target always goes in backward hemisphere
  G4bool incidentHasChanged = false;
  G4bool targetHasChanged = false;
  G4bool quasiElastic = false;
  G4FastVector<G4ReactionProduct,GHADLISTSIZE> vec;  // vec will contain the secondary particles
  G4int vecLen = 0;
  vec.Initialize(0);
    
  Cascade(vec, vecLen, originalIncident, currentParticle, targetParticle,
          incidentHasChanged, targetHasChanged, quasiElastic);
    
  CalculateMomenta(vec, vecLen, originalIncident, originalTarget,
                   modifiedOriginal, targetNucleus, currentParticle,
                   targetParticle, incidentHasChanged, targetHasChanged,
                   quasiElastic);
    
  SetUpChange(vec, vecLen, currentParticle, targetParticle, incidentHasChanged);
   
  if (isotopeProduction) DoIsotopeCounting(originalIncident, targetNucleus); 
  delete originalTarget;
  return &theParticleChange;
}

 
void G4LENeutronInelastic::SlowNeutron(const G4HadProjectile* originalIncident,
                                       G4ReactionProduct& modifiedOriginal,
                                       G4ReactionProduct& targetParticle,
                                       G4Nucleus& targetNucleus)
{        
  const G4double A = targetNucleus.GetA_asInt();    // atomic weight
  const G4double Z = targetNucleus.GetZ_asInt();    // atomic number
    
  G4double currentKinetic = modifiedOriginal.GetKineticEnergy()/MeV;
  G4double currentMass = modifiedOriginal.GetMass()/MeV;
  if( A < 1.5 )   // Hydrogen
    {
      //
      // very simple simulation of scattering angle and energy
      // nonrelativistic approximation with isotropic angular
      // distribution in the cms system 
      //
      G4double cost1, eka = 0.0;
      while (eka <= 0.0)
      {
        cost1 = -1.0 + 2.0*G4UniformRand();
        eka = 1.0 + 2.0*cost1*A + A*A;
      }
      G4double cost = std::min( 1.0, std::max( -1.0, (A*cost1+1.0)/std::sqrt(eka) ) );
      eka /= (1.0+A)*(1.0+A);
      G4double ek = currentKinetic*MeV/GeV;
      G4double amas = currentMass*MeV/GeV;
      ek *= eka;
      G4double en = ek + amas;
      G4double p = std::sqrt(std::abs(en*en-amas*amas));
      G4double sint = std::sqrt(std::abs(1.0-cost*cost));
      G4double phi = G4UniformRand()*twopi;
      G4double px = sint*std::sin(phi);
      G4double py = sint*std::cos(phi);
      G4double pz = cost;
      targetParticle.SetMomentum( px*GeV, py*GeV, pz*GeV );
      G4double pxO = originalIncident->Get4Momentum().x()/GeV;
      G4double pyO = originalIncident->Get4Momentum().y()/GeV;
      G4double pzO = originalIncident->Get4Momentum().z()/GeV;
      G4double ptO = pxO*pxO + pyO+pyO;
      if( ptO > 0.0 )
      {
        G4double pO = std::sqrt(pxO*pxO+pyO*pyO+pzO*pzO);
        cost = pzO/pO;
        sint = 0.5*(std::sqrt(std::abs((1.0-cost)*(1.0+cost)))+std::sqrt(ptO)/pO);
        G4double ph = pi/2.0;
        if( pyO < 0.0 )ph = ph*1.5;
        if( std::abs(pxO) > 0.000001 )ph = std::atan2(pyO,pxO);
        G4double cosp = std::cos(ph);
        G4double sinp = std::sin(ph);
        px = cost*cosp*px - sinp*py+sint*cosp*pz;
        py = cost*sinp*px + cosp*py+sint*sinp*pz;
        pz = -sint*px     + cost*pz;
      }
      else
      {
        if( pz < 0.0 )pz *= -1.0;
      }
      G4double pu = std::sqrt(px*px+py*py+pz*pz);
      modifiedOriginal.SetMomentum( targetParticle.GetMomentum() * (p/pu) );
      modifiedOriginal.SetKineticEnergy( ek*GeV );
      
      targetParticle.SetMomentum(
       originalIncident->Get4Momentum().vect() - modifiedOriginal.GetMomentum() );
      G4double pp = targetParticle.GetMomentum().mag();
      G4double tarmas = targetParticle.GetMass();
      targetParticle.SetTotalEnergy( std::sqrt( pp*pp + tarmas*tarmas ) );
      
      theParticleChange.SetEnergyChange( modifiedOriginal.GetKineticEnergy() );
      G4DynamicParticle *pd = new G4DynamicParticle;
      pd->SetDefinition( targetParticle.GetDefinition() );
      pd->SetMomentum( targetParticle.GetMomentum() );
      theParticleChange.AddSecondary( pd );
      return;
    }
    G4FastVector<G4ReactionProduct,4> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    G4double theAtomicMass = targetNucleus.AtomicMass( A, Z );
    G4double massVec[9];
    massVec[0] = targetNucleus.AtomicMass( A+1.0, Z     );
    massVec[1] = theAtomicMass;
    massVec[2] = 0.;
    if (Z > 1.0) 
        massVec[2] = targetNucleus.AtomicMass( A    , Z-1.0 );
    massVec[3] = 0.;
    if (Z > 1.0 && A > 1.0) 
        massVec[3] = targetNucleus.AtomicMass( A-1.0, Z-1.0 );
    massVec[4] = 0.;
    if (Z > 1.0 && A > 2.0 && A-2.0 > Z-1.0)
        massVec[4] = targetNucleus.AtomicMass( A-2.0, Z-1.0 );
    massVec[5] = 0.;
    if (Z > 2.0 && A > 3.0 && A-3.0 > Z-2.0)
        massVec[5] = targetNucleus.AtomicMass( A-3.0, Z-2.0 );
    massVec[6] = 0.;
    if (A > 1.0 && A-1.0 > Z) 
        massVec[6] = targetNucleus.AtomicMass( A-1.0, Z     );
    massVec[7] = massVec[3];
    massVec[8] = 0.;
    if (Z > 2.0 && A > 1.0)
        massVec[8] = targetNucleus.AtomicMass( A-1.0, Z-2.0 );
    
    theReactionDynamics.NuclearReaction( vec, vecLen, originalIncident,
                                         targetNucleus, theAtomicMass, massVec );
    
    theParticleChange.SetStatusChange( stopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );
    
    G4DynamicParticle * pd;
    for( G4int i=0; i<vecLen; ++i )
    {
      pd = new G4DynamicParticle();
      pd->SetDefinition( vec[i]->GetDefinition() );
      pd->SetMomentum( vec[i]->GetMomentum() );
      theParticleChange.AddSecondary( pd );
      delete vec[i];
    }
}


void G4LENeutronInelastic::Cascade(
   G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
   G4int& vecLen,
   const G4HadProjectile *originalIncident,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool &quasiElastic )
  {
    // derived from original FORTRAN code CASN by H. Fesefeldt (13-Sep-1987)
    //
    // Neutron undergoes interaction with nucleon within a nucleus.  Check if it is
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
    G4double centerofmassEnergy = std::sqrt( mOriginal*mOriginal +
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
    // npos = number of pi+, nneg = number of pi-, nzero = number of pi0
    G4int counter, nt=0, npos=0, nneg=0, nzero=0;
    const G4double c = 1.25;    
    const G4double b[] = { 0.35, 0.0 };
    if( first )      // compute normalization constants, this will only be Done once
    {
      first = false;
      G4int i;
      for( i=0; i<numMul; ++i )protmul[i] = 0.0;
      for( i=0; i<numSec; ++i )protnorm[i] = 0.0;
      counter = -1;
      for( npos=0; npos<numSec/3; ++npos )
      {
        for( nneg=std::max(0,npos-1); nneg<=(npos+1); ++nneg )
        {
          for( nzero=0; nzero<numSec/3; ++nzero )
          {
            if( ++counter < numMul )
            {
              nt = npos+nneg+nzero;
              if( nt > 0 )
              {
                protmul[counter] = Pmltpc(npos,nneg,nzero,nt,b[0],c) /
                  ( theReactionDynamics.Factorial(1-npos+nneg)*
                    theReactionDynamics.Factorial(1+npos-nneg) );
                protnorm[nt-1] += protmul[counter];
              }
            }
          }
        }
      }
      for( i=0; i<numMul; ++i )neutmul[i] = 0.0;
      for( i=0; i<numSec; ++i )neutnorm[i] = 0.0;
      counter = -1;
      for( npos=0; npos<(numSec/3); ++npos )
      {
        for( nneg=npos; nneg<=(npos+2); ++nneg )
        {
          for( nzero=0; nzero<numSec/3; ++nzero )
          {
            if( ++counter < numMul )
            {
              nt = npos+nneg+nzero;
              if( (nt>0) && (nt<=numSec) )
              {
                neutmul[counter] = Pmltpc(npos,nneg,nzero,nt,b[1],c) /
                  ( theReactionDynamics.Factorial(nneg-npos)*
                    theReactionDynamics.Factorial(2-nneg+npos) );
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
    
    const G4double expxu = 82.;      // upper bound for arg. of exp
    const G4double expxl = -expxu;        // lower bound for arg. of exp
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4int ieab = static_cast<G4int>(availableEnergy*5.0/GeV);
    const G4double supp[] = {0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98};
    G4double test, w0, wp, wt, wm;
    if( (availableEnergy < 2.0*GeV) && (G4UniformRand() >= supp[ieab]) )
    {
      // suppress high multiplicity events at low momentum
      // only one pion will be produced

      nneg = npos = nzero = 0;
      if( targetParticle.GetDefinition() == aNeutron )
      {
        test = std::exp( std::min( expxu, std::max( expxl, -(1.0+b[1])*(1.0+b[1])/(2.0*c*c) ) ) );
        w0 = test/2.0;
        wm = test;
        if( G4UniformRand() < w0/(w0+wm) )
          nzero = 1;
        else
          nneg = 1;
      } else { // target is a proton
        test = std::exp( std::min( expxu, std::max( expxl, -(1.0+b[0])*(1.0+b[0])/(2.0*c*c) ) ) );
        w0 = test;
        wp = test/2.0;        
        test = std::exp( std::min( expxu, std::max( expxl, -(-1.0+b[0])*(-1.0+b[0])/(2.0*c*c) ) ) );
        wm = test/2.0;
        wt = w0+wp+wm;
        wp += w0;
        G4double ran = G4UniformRand();
        if( ran < w0/wt )
          nzero = 1;
        else if( ran < wp/wt )
          npos = 1;
        else
          nneg = 1;
      }
    } else {  // (availableEnergy >= 2.0*GeV) || (random number < supp[ieab])
      G4double n, anpn;
      GetNormalizationConstant( availableEnergy, n, anpn );
      G4double ran = G4UniformRand();
      G4double dum, excs = 0.0;
      if( targetParticle.GetDefinition() == aProton )
      {
        counter = -1;
        for( npos=0; npos<numSec/3 && ran>=excs; ++npos )
        {
          for( nneg=std::max(0,npos-1); nneg<=(npos+1) && ran>=excs; ++nneg )
          {
            for( nzero=0; nzero<numSec/3 && ran>=excs; ++nzero )
            {
              if( ++counter < numMul )
              {
                nt = npos+nneg+nzero;
                if( nt > 0 )
                {
                  test = std::exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
                  if( std::fabs(dum) < 1.0 ) {
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
        npos--; nneg--; nzero--;
      } else { // target must be a neutron
        counter = -1;
        for( npos=0; npos<numSec/3 && ran>=excs; ++npos )
        {
          for( nneg=npos; nneg<=(npos+2) && ran>=excs; ++nneg )
          {
            for( nzero=0; nzero<numSec/3 && ran>=excs; ++nzero )
            {
              if( ++counter < numMul )
              {
                nt = npos+nneg+nzero;
                if( (nt>=1) && (nt<=numSec) )
                {
                  test = std::exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(nt*nt)/(n*n) ) ) );
                  dum = (pi/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
                  if( std::fabs(dum) < 1.0 ) {
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
        npos--; nneg--; nzero--;
      }
    }
    if( targetParticle.GetDefinition() == aProton )
    {
      switch( npos-nneg )
      {
       case 0:
         if( G4UniformRand() < 0.33 )
         {
           currentParticle.SetDefinitionAndUpdateE( aProton );
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
         currentParticle.SetDefinitionAndUpdateE( aProton );
         incidentHasChanged = true;
         break;
      }
    } else { // target must be a neutron
      switch( npos-nneg )
      {
       case -1:                       // changed from +1 by JLC, 7Jul97
         if( G4UniformRand() < 0.5 )
         {
           currentParticle.SetDefinitionAndUpdateE( aProton );
           incidentHasChanged = true;
         } else {
           targetParticle.SetDefinitionAndUpdateE( aProton );
           targetHasChanged = true;
         }
         break;
       case 0:
         break;
       default:
         currentParticle.SetDefinitionAndUpdateE( aProton );
         targetParticle.SetDefinitionAndUpdateE( aProton );
         incidentHasChanged = true;
         targetHasChanged = true;
         break;
      }
    }
    SetUpPions( npos, nneg, nzero, vec, vecLen );
// DEBUG -->    DumpFrames::DumpFrame(vec, vecLen);
    return;
  }

 /* end of file */
 
