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
// $Id: G4RPGInelastic.cc,v 1.2 2007/08/15 20:38:25 dennis Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#include "G4RPGInelastic.hh"
#include "Randomize.hh"
#include "G4HadReentrentException.hh"
#include "G4RPGStrangeProduction.hh"
#include "G4RPGTwoBody.hh"

 
G4double G4RPGInelastic::Pmltpc(G4int np, G4int nm, G4int nz, 
                                G4int n, G4double b, G4double c)
{
  const G4double expxu =  82.;           // upper bound for arg. of exp
  const G4double expxl = -expxu;         // lower bound for arg. of exp
  G4double npf = 0.0;
  G4double nmf = 0.0;
  G4double nzf = 0.0;
  G4int i;
  for( i=2; i<=np; i++ )npf += std::log((double)i);
  for( i=2; i<=nm; i++ )nmf += std::log((double)i);
  for( i=2; i<=nz; i++ )nzf += std::log((double)i);
  G4double r;
  r = std::min( expxu, std::max( expxl, -(np-nm+nz+b)*(np-nm+nz+b)/(2*c*c*n*n)-npf-nmf-nzf ) );
  return std::exp(r);
}


G4int G4RPGInelastic::Factorial( G4int n )
{
  G4int m = std::min(n,10);
  G4int result = 1;
  if( m <= 1 )return result;
  for( G4int i=2; i<=m; ++i )result *= i;
  return result;
}


G4bool G4RPGInelastic::MarkLeadingStrangeParticle(
     const G4ReactionProduct &currentParticle,
     const G4ReactionProduct &targetParticle,
     G4ReactionProduct &leadParticle )
{
  // The following was in GenerateXandPt and TwoCluster.
  // Add a parameter to the GenerateXandPt function telling it about the 
  // strange particle.
  //
  // Assumes that the original particle was a strange particle
  //
  G4bool lead = false;
  if( (currentParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
      (currentParticle.GetDefinition() != G4Proton::Proton()) &&
      (currentParticle.GetDefinition() != G4Neutron::Neutron()) )
  {
    lead = true;
    leadParticle = currentParticle;              //  set lead to the incident particle
  }
  else if( (targetParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
           (targetParticle.GetDefinition() != G4Proton::Proton()) &&
           (targetParticle.GetDefinition() != G4Neutron::Neutron()) )
  {
    lead = true;
    leadParticle = targetParticle;              //   set lead to the target particle
  }
  return lead;
}

 
 void G4RPGInelastic::SetUpPions(const G4int np, const G4int nm, 
                                 const G4int nz,
                             G4FastVector<G4ReactionProduct,256> &vec,
                                 G4int &vecLen)
 {
   if( np+nm+nz == 0 )return;
   G4int i;
   G4ReactionProduct *p;
   for( i=0; i<np; ++i )
   {
     p = new G4ReactionProduct;
     p->SetDefinition( G4PionPlus::PionPlus() );
     (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
     vec.SetElement( vecLen++, p );
   }
   for( i=np; i<np+nm; ++i )
   {
     p = new G4ReactionProduct;
     p->SetDefinition( G4PionMinus::PionMinus() );
     (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
     vec.SetElement( vecLen++, p );
   }
   for( i=np+nm; i<np+nm+nz; ++i )
   {
     p = new G4ReactionProduct;
     p->SetDefinition( G4PionZero::PionZero() );
     (G4UniformRand() < 0.5) ? p->SetSide( -1 ) : p->SetSide( 1 );
     vec.SetElement( vecLen++, p );
   }
 }


 void G4RPGInelastic::GetNormalizationConstant(
   const G4double energy,  // MeV,  <0 means annihilation channels
   G4double &n,
   G4double &anpn )
  {
    const G4double expxu =  82.;          // upper bound for arg. of exp
    const G4double expxl = -expxu;        // lower bound for arg. of exp
    const G4int numSec = 60;
    //
    // the only difference between the calculation for annihilation channels
    // and normal is the starting value, iBegin, for the loop below
    //
    G4int iBegin = 1;
    G4double en = energy;
    if( energy < 0.0 )
    {
      iBegin = 2;
      en *= -1.0;
    }
    //
    // number of total particles vs. centre of mass Energy - 2*proton mass
    //
    G4double aleab = std::log(en/GeV);
    n = 3.62567 + aleab*(0.665843 + aleab*(0.336514 + aleab*(0.117712 + 0.0136912*aleab)));
    n -= 2.0;
    //
    // normalization constant for kno-distribution
    //
    anpn = 0.0;
    G4double test, temp;
    for( G4int i=iBegin; i<=numSec; ++i )
    {
      temp = pi*i/(2.0*n*n);
      test = std::exp( std::min( expxu, std::max( expxl, -(pi/4.0)*(i*i)/(n*n) ) ) );
      if( temp < 1.0 )
      {
        if( test >= 1.0e-10 )anpn += temp*test;
      }
      else
        anpn += temp*test;
    }
  }
 
void G4RPGInelastic::CalculateMomenta(
  G4FastVector<G4ReactionProduct,256> &vec,
  G4int &vecLen,
  const G4HadProjectile *originalIncident,
  const G4DynamicParticle *originalTarget,
  G4ReactionProduct &modifiedOriginal,
  G4Nucleus &targetNucleus,
  G4ReactionProduct &currentParticle,
  G4ReactionProduct &targetParticle,
  G4bool &incidentHasChanged,
  G4bool &targetHasChanged,
  G4bool quasiElastic )
{
  cache = 0;
  what = originalIncident->Get4Momentum().vect();

  G4ReactionProduct leadingStrangeParticle;

  strangeProduction.ReactionStage(originalIncident, modifiedOriginal,
                                  incidentHasChanged, originalTarget,
                                  targetParticle, targetHasChanged,
                                  targetNucleus, currentParticle,
                                  vec, vecLen,
				  false, leadingStrangeParticle);

  if( quasiElastic )
  {
    twoBody.ReactionStage(originalIncident, modifiedOriginal,
                          incidentHasChanged, originalTarget,
                          targetParticle, targetHasChanged,
                          targetNucleus, currentParticle,
                          vec, vecLen,
                          false, leadingStrangeParticle);
    return;
  }

  G4bool leadFlag = MarkLeadingStrangeParticle(currentParticle,
                                               targetParticle,
                                               leadingStrangeParticle );
  //
  // Note: the number of secondaries can be reduced in GenerateXandPt 
  // and TwoCluster
  //
  G4bool finishedGenXPt = false;
  G4bool annihilation = false;
  if( originalIncident->GetDefinition()->GetPDGEncoding() < 0 &&
      currentParticle.GetMass() == 0.0 && targetParticle.GetMass() == 0.0 )
  {
    // original was an anti-particle and annihilation has taken place
    annihilation = true;
    G4double ekcor = 1.0;
    G4double ek = originalIncident->GetKineticEnergy();
    G4double ekOrg = ek;
      
    const G4double tarmas = originalTarget->GetDefinition()->GetPDGMass();
    if( ek > 1.0*GeV )ekcor = 1./(ek/GeV);
    const G4double atomicWeight = targetNucleus.GetN();
    ek = 2*tarmas + ek*(1.+ekcor/atomicWeight);
    G4double tkin = targetNucleus.Cinema(ek);
    ek += tkin;
    ekOrg += tkin;
    //      modifiedOriginal.SetKineticEnergy( ekOrg );
    //
    // evaporation --  re-calculate black track energies
    //                 this was Done already just before the cascade
    //
    tkin = targetNucleus.AnnihilationEvaporationEffects(ek, ekOrg);
    ekOrg -= tkin;
    ekOrg = std::max( 0.0001*GeV, ekOrg );
    modifiedOriginal.SetKineticEnergy( ekOrg );
    G4double amas = originalIncident->GetDefinition()->GetPDGMass();
    G4double et = ekOrg + amas;
    G4double p = std::sqrt( std::abs(et*et-amas*amas) );
    G4double pp = modifiedOriginal.GetMomentum().mag();
    if( pp > 0.0 )
    {
      G4ThreeVector momentum = modifiedOriginal.GetMomentum();
      modifiedOriginal.SetMomentum( momentum * (p/pp) );
    }
    if( ekOrg <= 0.0001 )
    {
      modifiedOriginal.SetKineticEnergy( 0.0 );
      modifiedOriginal.SetMomentum( 0.0, 0.0, 0.0 );
    }
  }

  // twsup gives percentage of time two-cluster model is called

  const G4double twsup[] = { 1.0, 0.7, 0.5, 0.3, 0.2, 0.1 };
  G4double rand1 = G4UniformRand();
  G4double rand2 = G4UniformRand();

  // Cache current, target, and secondaries
  G4ReactionProduct saveCurrent = currentParticle;
  G4ReactionProduct saveTarget = targetParticle;
  std::vector<G4ReactionProduct> savevec;
  for (G4int i = 0; i < vecLen; i++) savevec.push_back(*vec[i]);

  if( annihilation || (vecLen >= 6) ||
      (modifiedOriginal.GetKineticEnergy()/GeV >= 1.0) &&
      (((originalIncident->GetDefinition() == G4KaonPlus::KaonPlus() ||
         originalIncident->GetDefinition() == G4KaonMinus::KaonMinus() ||
         originalIncident->GetDefinition() == G4KaonZeroLong::KaonZeroLong() ||
         originalIncident->GetDefinition() == G4KaonZeroShort::KaonZeroShort()) &&
        rand1 < 0.5) || rand2 > twsup[vecLen]) )

    finishedGenXPt =
      fragmentation.ReactionStage(originalIncident, modifiedOriginal,
                                  incidentHasChanged, originalTarget,
                                  targetParticle, targetHasChanged,
                                  targetNucleus, currentParticle,
                                  vec, vecLen,
                                  leadFlag, leadingStrangeParticle);

  if( finishedGenXPt )
  {
    Rotate(vec, vecLen);
    return;
  }

  G4bool finishedTwoClu = false;
  if( modifiedOriginal.GetTotalMomentum()/MeV < 1.0 )
  {
    for(G4int i=0; i<vecLen; i++) delete vec[i];
    vecLen = 0;
  }
  else
  {
    // Occaisionally, GenerateXandPt will fail in the annihilation channel.
    // Restore current, target and secondaries to pre-GenerateXandPt state
    // before trying annihilation in TwoCluster

    if (!finishedGenXPt && annihilation) {
      currentParticle = saveCurrent;
      targetParticle = saveTarget;
      for (G4int i = 0; i < vecLen; i++) delete vec[i];
      vecLen = 0;
      vec.Initialize( 0 );
      for (G4int i = 0; i < G4int(savevec.size()); i++) {
        G4ReactionProduct* p = new G4ReactionProduct;
        *p = savevec[i];
        vec.SetElement( vecLen++, p );
      }
    }

    pionSuppression.ReactionStage(originalIncident, modifiedOriginal,
                                  incidentHasChanged, originalTarget,
                                  targetParticle, targetHasChanged,
                                  targetNucleus, currentParticle,
                                  vec, vecLen,
                                  false, leadingStrangeParticle);

    try
    {
      finishedTwoClu = 
        twoCluster.ReactionStage(originalIncident, modifiedOriginal,
                                 incidentHasChanged, originalTarget,
                                 targetParticle, targetHasChanged,
                                 targetNucleus, currentParticle,
                                 vec, vecLen,
                                 leadFlag, leadingStrangeParticle);
    }
     catch(G4HadReentrentException aC)
    {
       aC.Report(G4cout);
       throw G4HadReentrentException(__FILE__, __LINE__, "Failing to calculate momenta");
    }
  }

  if( finishedTwoClu )
  {
    Rotate(vec, vecLen);
    return;
  }

  twoBody.ReactionStage(originalIncident, modifiedOriginal,
                        incidentHasChanged, originalTarget,
                        targetParticle, targetHasChanged,
                        targetNucleus, currentParticle,
                        vec, vecLen,
                        false, leadingStrangeParticle);
}

 
 void G4RPGInelastic::
 Rotate(G4FastVector<G4ReactionProduct,256> &vec, G4int &vecLen)
 {
   G4double rotation = 2.*pi*G4UniformRand();
   cache = rotation;
   G4int i;
   for( i=0; i<vecLen; ++i )
   {
     G4ThreeVector momentum = vec[i]->GetMomentum();
     momentum = momentum.rotate(rotation, what);
     vec[i]->SetMomentum(momentum);
   }
 }      

void 
G4RPGInelastic::SetUpChange(G4FastVector<G4ReactionProduct,256> &vec,
                            G4int &vecLen,
                            G4ReactionProduct &currentParticle,
                            G4ReactionProduct &targetParticle,
                            G4bool &incidentHasChanged )
{
  theParticleChange.Clear();
  G4ParticleDefinition *aKaonZL = G4KaonZeroLong::KaonZeroLong();
  G4ParticleDefinition *aKaonZS = G4KaonZeroShort::KaonZeroShort();
  G4int i;
  if( currentParticle.GetDefinition() == aKaonZL )
  {
    if( G4UniformRand() <= 0.5 )
    {
      currentParticle.SetDefinition( aKaonZS );
      incidentHasChanged = true;
    }
  }
  else if( currentParticle.GetDefinition() == aKaonZS )
  {
    if( G4UniformRand() > 0.5 )
    {
      currentParticle.SetDefinition( aKaonZL );
      incidentHasChanged = true;
    }
  }

  if( targetParticle.GetDefinition() == aKaonZL )
  {
    if( G4UniformRand() <= 0.5 )targetParticle.SetDefinition( aKaonZS );
  }
  else if( targetParticle.GetDefinition() == aKaonZS )
  {
    if( G4UniformRand() > 0.5 )targetParticle.SetDefinition( aKaonZL );
  }
  for( i=0; i<vecLen; ++i )
  {
    if( vec[i]->GetDefinition() == aKaonZL )
    {
      if( G4UniformRand() <= 0.5 )vec[i]->SetDefinition( aKaonZS );
    }
    else if( vec[i]->GetDefinition() == aKaonZS )
    {
      if( G4UniformRand() > 0.5 )vec[i]->SetDefinition( aKaonZL );
    }
  }

  if( incidentHasChanged )
  {
    G4DynamicParticle* p0 = new G4DynamicParticle;
    p0->SetDefinition( currentParticle.GetDefinition() );
    p0->SetMomentum( currentParticle.GetMomentum() );
    theParticleChange.AddSecondary( p0 );
    theParticleChange.SetStatusChange( stopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );
  }
  else
  {
    G4double p = currentParticle.GetMomentum().mag()/MeV;
    G4ThreeVector m = currentParticle.GetMomentum();
    if( p > DBL_MIN )
      theParticleChange.SetMomentumChange( m.x()/p, m.y()/p, m.z()/p );
    else
      theParticleChange.SetMomentumChange( 0.0, 0.0, 1.0 );

    G4double aE = currentParticle.GetKineticEnergy();
    if (std::fabs(aE)<.1*eV) aE=.1*eV;
    theParticleChange.SetEnergyChange( aE );
  }

  if( targetParticle.GetMass() > 0.0 )  // Tgt particle can be eliminated in TwoBody
  {
    G4ThreeVector momentum = targetParticle.GetMomentum();
    momentum = momentum.rotate(cache, what);
    G4double targKE = targetParticle.GetKineticEnergy();
    G4ThreeVector dir(0.0, 0.0, 1.0);
    if (targKE < DBL_MIN)
      targKE = DBL_MIN;
    else
      dir = momentum/momentum.mag();

    G4DynamicParticle* p1 = 
      new G4DynamicParticle(targetParticle.GetDefinition(), dir, targKE);

    theParticleChange.AddSecondary( p1 );
  }

  G4DynamicParticle* p;
  for( i=0; i<vecLen; ++i )
  {
    G4double secKE = vec[i]->GetKineticEnergy();
    G4ThreeVector momentum = vec[i]->GetMomentum();
    G4ThreeVector dir(0.0, 0.0, 1.0);
    if (secKE < DBL_MIN)
      secKE = DBL_MIN;
    else
      dir = momentum/momentum.mag();

    p = new G4DynamicParticle(vec[i]->GetDefinition(), dir, secKE);
    theParticleChange.AddSecondary( p );
    delete vec[i];
  }
}
 
/* end of file */
