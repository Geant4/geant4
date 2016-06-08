// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
 // Hadronic Process: Inelastic Interaction 
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 27-Mar-1997
 // J.P. Wellisch: 23-Apr-97: G4Exception removed
 // J.P. Wellisch: 24-Apr-97: correction for SetUpPions
 // Modified by J.L. Chuma, 30-Apr-97: added originalTarget to CalculateMomenta
 //                                    since TwoBody needed to reset the target particle
 // J.L. Chuma, 20-Jun-97: Modified CalculateMomenta to correct the decision process
 //                        for whether to use GenerateXandPt or TwoCluster
 // J.L. Chuma, 06-Aug-97: added original incident particle, before Fermi motion and
 //                        evaporation effects are included, needed for calculating
 //                        self absorption and corrections for single particle spectra
 // HPW removed misunderstanding of LocalEnergyDeposit, 11.04.98.
 
#include "G4InelasticInteraction.hh"
#include "Randomize.hh"
 
 G4double
  G4InelasticInteraction::Pmltpc(      // used in Cascade functions
   G4int np, G4int nm, G4int nz, G4int n, G4double b, G4double c )
  {
    const G4double expxu =  82.;           // upper bound for arg. of exp
    const G4double expxl = -expxu;         // lower bound for arg. of exp
    G4double npf = 0.0;
    G4double nmf = 0.0;
    G4double nzf = 0.0;
    G4int i;
    for( i=2; i<=np; i++ )npf += log((double)i);
    for( i=2; i<=nm; i++ )nmf += log((double)i);
    for( i=2; i<=nz; i++ )nzf += log((double)i);
    G4double r;
    r = G4std::min( expxu, G4std::max( expxl, -(np-nm+nz+b)*(np-nm+nz+b)/(2*c*c*n*n)-npf-nmf-nzf ) );
    return exp(r);
  }

 G4bool
  G4InelasticInteraction::MarkLeadingStrangeParticle(
   const G4ReactionProduct &currentParticle,
   const G4ReactionProduct &targetParticle,
   G4ReactionProduct &leadParticle )
  {
    // the following was in GenerateXandPt and TwoCluster
    // add a parameter to the GenerateXandPt function telling it about the strange particle
    //
    // assumes that the original particle was a strange particle
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
 
 void
  G4InelasticInteraction::SetUpPions(
   const G4int np,
   const G4int nm,
   const G4int nz,
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int &vecLen )
  {
    if( np+nm+nz == 0 )return;
    G4int i;
    G4ReactionProduct *p = new G4ReactionProduct [np+nm+nz];
    for( i=0; i<np; ++i )
    {
      p[i].SetDefinition( G4PionPlus::PionPlus() );
      (G4UniformRand() < 0.5) ? p[i].SetSide( -1 ) : p[i].SetSide( 1 );
      vec.SetElement( vecLen++, &p[i] );
    }
    for( i=np; i<np+nm; ++i )
    {
      p[i].SetDefinition( G4PionMinus::PionMinus() );
      (G4UniformRand() < 0.5) ? p[i].SetSide( -1 ) : p[i].SetSide( 1 );
      vec.SetElement( vecLen++, &p[i] );
    }
    for( i=np+nm; i<np+nm+nz; ++i )
    {
      p[i].SetDefinition( G4PionZero::PionZero() );
      (G4UniformRand() < 0.5) ? p[i].SetSide( -1 ) : p[i].SetSide( 1 );
      vec.SetElement( vecLen++, &p[i] );
    }
  }
 
 void
  G4InelasticInteraction::GetNormalizationConstant(
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
    G4double aleab = log(en/GeV);
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
      test = exp( G4std::min( expxu, G4std::max( expxl, -(pi/4.0)*(i*i)/(n*n) ) ) );
      if( temp < 1.0 )
      {
        if( test >= 1.0e-10 )anpn += temp*test;
      }
      else
        anpn += temp*test;
    }
  }
 
 void
  G4InelasticInteraction::CalculateMomenta(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int &vecLen,
   const G4DynamicParticle *originalIncident,   // the original incident particle
   const G4DynamicParticle *originalTarget,
   G4ReactionProduct &modifiedOriginal,   // Fermi motion and evap. effects included
   G4Nucleus &targetNucleus,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged,
   G4bool &targetHasChanged,
   G4bool quasiElastic )
  {
    theReactionDynamics.ProduceStrangeParticlePairs( vec, vecLen,
                                                     modifiedOriginal, originalTarget,
                                                     currentParticle, targetParticle,
                                                     incidentHasChanged, targetHasChanged );
    if( quasiElastic )
    {
      theReactionDynamics.TwoBody( vec, vecLen,
                                   modifiedOriginal, originalTarget,
                                   currentParticle, targetParticle,
                                   targetNucleus, targetHasChanged );
      return;
    }
    G4ReactionProduct leadingStrangeParticle;
    G4bool leadFlag = MarkLeadingStrangeParticle( currentParticle,
                                                  targetParticle,
                                                  leadingStrangeParticle );
    //
    // Note: the number of secondaries can be reduced in GenerateXandPt and TwoCluster
    //
    G4bool finishedGenXPt = false;
    G4bool annihilation = false;
    if( originalIncident->GetDefinition()->GetPDGEncoding() < 0 &&
        currentParticle.GetMass() == 0.0 && targetParticle.GetMass() == 0.0 )
    {
      // original was an anti-particle and annihilation has taken place
      annihilation = true;
      G4double ekcor = 1.0;
      G4double ek = originalIncident->GetKineticEnergy()/GeV;
      const G4double tarmas = originalTarget->GetDefinition()->GetPDGMass()/GeV;
      if( ek > 1.0 )ekcor = 1./ek;
      const G4double atomicWeight = targetNucleus.GetN();
      ek = 2*tarmas + ek*(1.+ekcor/atomicWeight);
      modifiedOriginal.SetKineticEnergy( ek*GeV );
      //
      // evaporation --  re-calculate black track energies
      //                 this was Done already just before the cascade
      //
      G4double tkin = targetNucleus.EvaporationEffects( ek*GeV )/GeV;
      ek -= tkin;
      ek = G4std::max( 0.0001, ek );
      modifiedOriginal.SetKineticEnergy( ek*GeV );
      G4double amas = originalIncident->GetDefinition()->GetPDGMass()/GeV;
      G4double et = ek + amas;
      G4double p = sqrt( abs(et*et-amas*amas) );
      G4double pp = modifiedOriginal.GetMomentum().mag()/GeV;
      if( pp > 0.0 )
      {
        G4ThreeVector momentum = modifiedOriginal.GetMomentum();
        modifiedOriginal.SetMomentum( momentum * (p/pp) );
      }
      if( ek <= 0.0001 )
      {
        modifiedOriginal.SetKineticEnergy( 0.0 );
        modifiedOriginal.SetMomentum( 0.0, 0.0, 0.0 );
      }
    }
    const G4double twsup[] = { 1.0, 0.7, 0.5, 0.3, 0.2, 0.1 };
    G4double rand1 = G4UniformRand();
    G4double rand2 = G4UniformRand();
    if( annihilation || (vecLen >= 6) ||
        (modifiedOriginal.GetKineticEnergy()/GeV >= 1.0) &&
        (((originalIncident->GetDefinition() == G4KaonPlus::KaonPlus() ||
           originalIncident->GetDefinition() == G4KaonMinus::KaonMinus() ||
           originalIncident->GetDefinition() == G4KaonZeroLong::KaonZeroLong() ||
           originalIncident->GetDefinition() == G4KaonZeroShort::KaonZeroShort()) &&
          rand1 < 0.5) || rand2 > twsup[vecLen]) )
      finishedGenXPt =
        theReactionDynamics.GenerateXandPt( vec, vecLen,
                                            modifiedOriginal, originalIncident,
                                            currentParticle, targetParticle,
                                            targetNucleus, incidentHasChanged,
                                            targetHasChanged, leadFlag,
                                            leadingStrangeParticle );
    if( finishedGenXPt )return;
    G4bool finishedTwoClu = false;
    if( modifiedOriginal.GetTotalMomentum()/MeV < 1.0 )vecLen = 0;
    else
    {
      theReactionDynamics.SuppressChargedPions( vec, vecLen,
                                                modifiedOriginal, currentParticle,
                                                targetParticle, targetNucleus,
                                                incidentHasChanged, targetHasChanged );
      finishedTwoClu = theReactionDynamics.TwoCluster( vec, vecLen,
                                                       modifiedOriginal, originalIncident,
                                                       currentParticle, targetParticle,
                                                       targetNucleus, incidentHasChanged,
                                                       targetHasChanged, leadFlag,
                                                       leadingStrangeParticle );
    }
    if( finishedTwoClu )return;
    //
    // PNBlackTrackEnergy is the kinetic energy available for
    //   proton/neutron black track particles [was enp(1) in fortran code]
    // DTABlackTrackEnergy is the kinetic energy available for
    //   deuteron/triton/alpha particles      [was enp(3) in fortran code]
    //const G4double pnCutOff = 0.1;
    //const G4double dtaCutOff = 0.1;
    //if( (targetNucleus.GetN() >= 1.5)
    //    && !(incidentHasChanged || targetHasChanged)
    //    && (targetNucleus.GetPNBlackTrackEnergy()/MeV <= pnCutOff)
    //    && (targetNucleus.GetDTABlackTrackEnergy()/MeV <= dtaCutOff) )
    //{
    // the atomic weight of the target nucleus is >= 1.5            AND
    //   neither the incident nor the target particles have changed  AND
    //     there is no kinetic energy available for either proton/neutron
    //     or for deuteron/triton/alpha black track particles
    // For diffraction scattering on heavy nuclei use elastic routines instead
    //G4cerr << "*** Error in G4InelasticInteraction::CalculateMomenta" << G4endl;
    //G4cerr << "*** the elastic scattering would be better here ***" <<G4endl;
    //}
    theReactionDynamics.TwoBody( vec, vecLen,
                                 modifiedOriginal, originalTarget,
                                 currentParticle, targetParticle,
                                 targetNucleus, targetHasChanged );
  }
 
 void
  G4InelasticInteraction::SetUpChange(
   G4FastVector<G4ReactionProduct,128> &vec,
   G4int &vecLen,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4bool &incidentHasChanged )
  {
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
      theParticleChange.SetNumberOfSecondaries( vecLen+2 );
      G4DynamicParticle* p0 = new G4DynamicParticle;
      p0->SetDefinition( currentParticle.GetDefinition() );
      p0->SetMomentum( currentParticle.GetMomentum() );
      theParticleChange.AddSecondary( p0 );
      theParticleChange.SetStatusChange( fStopAndKill );
      theParticleChange.SetEnergyChange( 0.0 );
    }
    else
    {
      theParticleChange.SetNumberOfSecondaries( vecLen+1 );
      G4double p = currentParticle.GetMomentum().mag()/MeV;
      G4ThreeVector m = currentParticle.GetMomentum();
      if( p > DBL_MIN )
        theParticleChange.SetMomentumChange( m.x()/p, m.y()/p, m.z()/p );
      else
        theParticleChange.SetMomentumChange( 0.0, 0.0, 0.0 );

      theParticleChange.SetEnergyChange( currentParticle.GetKineticEnergy() );
    }
    if( targetParticle.GetMass() > 0.0 )  // targetParticle can be eliminated in TwoBody
    {
      G4DynamicParticle *p1 = new G4DynamicParticle;
      p1->SetDefinition( targetParticle.GetDefinition() );
      p1->SetMomentum( targetParticle.GetMomentum() );
      theParticleChange.AddSecondary( p1 );
    }
    G4DynamicParticle *p;
    for( i=0; i<vecLen; ++i )
    {
      p = new G4DynamicParticle();
      p->SetDefinition( vec[i]->GetDefinition() );
      p->SetMomentum( vec[i]->GetMomentum() );
      theParticleChange.AddSecondary( p );
    }
  }
 
  /* end of file */
 
