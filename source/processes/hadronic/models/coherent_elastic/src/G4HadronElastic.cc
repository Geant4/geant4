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
// $Id: G4HadronElastic.cc,v 1.68 2010-11-19 18:50:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Physics model class G4HadronElastic (derived from G4LElastic)
//
//
// G4 Model: Low-energy Elastic scattering with 4-momentum balance
// F.W. Jones, TRIUMF, 04-JUN-96
// Uses  G4ElasticHadrNucleusHE and G4VQCrossSection
//
//
// 25-JUN-98 FWJ: replaced missing Initialize for ParticleChange.
// 09-Set-05 V.Ivanchenko HARP version of the model: fix scattering
//           on hydrogen, use relativistic Lorentz transformation
// 24-Nov-05 V.Ivanchenko sample cost in center of mass reference system
// 03-Dec-05 V.Ivanchenko add protection to initial momentum 20 MeV/c in
//           center of mass system (before it was in lab system)
//           below model is not valid
// 14-Dec-05 V.Ivanchenko change protection to cos(theta) < -1 and
//           rename the class
// 13-Apr-06 V.Ivanchenko move to coherent_elastic subdirectory; remove
//           charge exchange; remove limitation on incident momentum;
//           add s-wave regim below some momentum        
// 24-Apr-06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07-Jun-06 V.Ivanchenko fix problem of rotation
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy, below which S-wave is sampled
// 02-Aug-06 V.Ivanchenko introduce energy cut on the aria of S-wave for pions
// 24-Aug-06 V.Ivanchenko switch on G4ElasticHadrNucleusHE
// 31-Aug-06 V.Ivanchenko do not sample sacttering for particles with kinetic 
//                        energy below 10 keV
// 16-Nov-06 V.Ivanchenko Simplify logic of choosing of the model for sampling
// 30-Mar-07 V.Ivanchenko lowEnergyLimitQ=0, lowEnergyLimitHE = 1.0*GeV,
//                        lowestEnergyLimit= 0
// 04-May-07 V.Ivanchenko do not use HE model for hydrogen target to avoid NaN;
//                        use QElastic for p, n incident for any energy for 
//                        p and He targets only  
// 11-May-07 V.Ivanchenko remove unused method Defs1
// 13.01.10: M.Kosov: Use G4Q(Pr/Neut)ElasticCS instead of G4QElasticCS
//

#include "G4HadronElastic.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"
#include "G4VQCrossSection.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

G4VQCrossSection* G4HadronElastic::pCManager = 0;
G4VQCrossSection* G4HadronElastic::nCManager = 0;

G4HadronElastic::G4HadronElastic(G4ElasticHadrNucleusHE* HModel) 
  : G4HadronicInteraction("G4HadronElastic"), hElastic(HModel)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  verboseLevel= 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ  = 0.0*GeV;  
  lowEnergyLimitHE = 1.0*GeV; 
  lowestEnergyLimit= 1.e-6*eV;  
  plabLowLimit     = 20.0*MeV;

  if(!pCManager)
  {
    pCManager = G4QProtonElasticCrossSection::GetPointer();
    nCManager = G4QNeutronElasticCrossSection::GetPointer();
  }
  if(!hElastic) hElastic = new G4ElasticHadrNucleusHE();

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();
  thePionPlus = G4PionPlus::PionPlus();
  thePionMinus= G4PionMinus::PionMinus();

}

G4HadronElastic::~G4HadronElastic()
{
  delete hElastic;
}

G4VQCrossSection* G4HadronElastic::GetCS()
{
  return pCManager;
  //if     (PDG==2212) return pCManager;
  //else if(PDG==2112) return nCManager;
  //return 0;
}

G4ElasticHadrNucleusHE* G4HadronElastic::GetHElastic()
{
  return hElastic;
}

G4HadFinalState* G4HadronElastic::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin <= lowestEnergyLimit) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  G4double plab = aParticle->GetTotalMomentum();
  if (verboseLevel >1) {
    G4cout << "G4HadronElastic::DoIt: Incident particle plab=" 
	   << plab/GeV << " GeV/c " 
	   << " ekin(MeV) = " << ekin/MeV << "  " 
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;
  }
  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();

  G4int N = A - Z;
  G4int projPDG = theParticle->GetPDGEncoding();
  if (verboseLevel>1) {
    G4cout << "G4HadronElastic for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z 
	   << " A= " << A << " N= " << N 
	   << G4endl;
  }
  G4ParticleDefinition * theDef = 0;

  if(Z == 1 && A == 1)       theDef = theProton;
  else if (Z == 1 && A == 2) theDef = theDeuteron;
  else if (Z == 1 && A == 3) theDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theDef = G4He3::He3();
  else if (Z == 2 && A == 4) theDef = theAlpha;
  else theDef = G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.0);
 
  G4double m2 = theDef->GetPDGMass();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,m2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();
  G4double tmax = 4.0*ptot*ptot;
  G4double t = 0.0;

  // Choose generator
  G4ElasticGenerator gtype = fLElastic;

  // Q-elastic for p,n scattering on H and He
  if (theParticle == theProton || theParticle == theNeutron) {
    //     && Z <= 2 && ekin >= lowEnergyLimitQ)  
    gtype = fQElastic;

  } else {
    // S-wave for very low energy
    if(plab < plabLowLimit) gtype = fSWave;
    // HE-elastic for energetic projectile mesons
    else if(ekin >= lowEnergyLimitHE && theParticle->GetBaryonNumber() == 0) 
      { gtype = fHElastic; }
  }

  //
  // Sample t
  //
  if(gtype == fQElastic) {
    if (verboseLevel >1) {
      G4cout << "G4HadronElastic: Z= " << Z << " N= "  << N << " pdg= " <<  projPDG
	       << " mom(GeV)= " << plab/GeV<<", pC="<<pCManager<<", nC="<<nCManager<<G4endl; 
    }
    if(Z == 1 && N == 2) N = 1;
    else if(Z == 2 && N == 1) N = 2;
    G4double cs = 0.;
    if     (projPDG==2212) cs = pCManager->GetCrossSection(false,plab,Z,N,projPDG);
    else if(projPDG==2112) cs = nCManager->GetCrossSection(false,plab,Z,N,projPDG);

    // check if cross section is reasonable
    if(cs > 0.0)
    {
      if     (projPDG==2212) t = pCManager->GetExchangeT(Z,N,projPDG);
      else if(projPDG==2112) t = nCManager->GetExchangeT(Z,N,projPDG);
    }
    else if(plab > plabLowLimit) gtype = fLElastic;
    else gtype = fSWave;
  }

  if(gtype == fLElastic) {
    G4double g2 = GeV*GeV; 
    t = g2*SampleT(tmax/g2,m1,m2, A);
  }

  // use mean atomic number
  if(gtype == fHElastic) {
    t = hElastic->SampleT(theParticle,plab, Z, A);
  }

  if(gtype == fSWave) t = G4UniformRand()*tmax;

  if(verboseLevel>1) {
    G4cout <<"type= " << gtype <<" t= " << t << " tmax= " << tmax 
	   << " ptot= " << ptot << G4endl;
  }
  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  // problem in sampling
  if(cost > 1.0 || cost < -1.0) {
    if(verboseLevel > 0) {
      G4cout << "G4HadronElastic:WARNING: Z= " << Z << " N= " 
	     << N << " " << aParticle->GetDefinition()->GetParticleName()
	     << " mom(GeV)= " << plab/GeV 
	     << " the model type " << gtype;
      if(gtype ==  fQElastic) G4cout << " CHIPS ";
      else if(gtype ==  fLElastic) G4cout << " LElastic ";
      else if(gtype ==  fHElastic) G4cout << " HElastic ";
      G4cout << " cost= " << cost 
	     << G4endl; 
    }
    cost = 1.0;
    sint = 0.0;

    // normal situation
  } else  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }    
  if (verboseLevel>1) {
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) {
    G4cout << "Scattered: "
	   << nlv1<<" m= " << m1 << " ekin(MeV)= " << eFinal 
	   << " Proj: 4-mom " << lv1 
	   <<G4endl;
  }
  if(eFinal <= lowestEnergyLimit) {
    if(eFinal < 0.0 && verboseLevel > 0) {
      G4cout << "G4HadronElastic WARNING ekin= " << eFinal
	     << " after scattering of " 
	     << aParticle->GetDefinition()->GetParticleName()
	     << " p(GeV/c)= " << plab
	     << " on " << theDef->GetParticleName()
	     << G4endl;
    }
    theParticleChange.SetEnergyChange(0.0);
    nlv1 = G4LorentzVector(0.0,0.0,0.0,m1);

  } else {
    theParticleChange.SetMomentumChange(nlv1.vect().unit());
    theParticleChange.SetEnergyChange(eFinal);
  }  

  G4LorentzVector nlv0 = lv - nlv1;
  G4double erec =  nlv0.e() - m2;
  if (verboseLevel > 1) {
    G4cout << "Recoil: "
	   << nlv0<<" m= " << m2 << " ekin(MeV)= " << erec 
	   <<G4endl;
  }
  if(erec > lowEnergyRecoilLimit) {
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, nlv0);
    theParticleChange.AddSecondary(aSec);
  } else {
    if(erec < 0.0) erec = 0.0;
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}

G4double 
G4HadronElastic::SampleT(G4double tmax, G4double, G4double, G4double atno2)
{
  // G4cout << "Entering elastic scattering 2"<<G4endl;
  // Compute the direction of elastic scattering.
  // It is planned to replace this code with a method based on
  // parameterized functions and a Monte Carlo method to invert the CDF.

  //  G4double ran = G4UniformRand();
  G4double aa, bb, cc, dd, rr;
  if (atno2 <= 62.) {
    aa = std::pow(atno2, 1.63);
    bb = 14.5*std::pow(atno2, 0.66);
    cc = 1.4*std::pow(atno2, 0.33);
    dd = 10.;
  } else {
    aa = std::pow(atno2, 1.33);
    bb = 60.*std::pow(atno2, 0.33);
    cc = 0.4*std::pow(atno2, 0.40);
    dd = 10.;
  }
  aa = aa/bb;
  cc = cc/dd;
  G4double ran, t1, t2;
  do {
    ran = G4UniformRand();
    t1 = -std::log(ran)/bb;
    t2 = -std::log(ran)/dd;
  } while(t1 > tmax || t2 > tmax);

  rr = (aa + cc)*ran;

  if (verboseLevel > 1) {
    G4cout << "DoIt: aa,bb,cc,dd,rr" << G4endl;
    G4cout << aa << " " << bb << " " << cc << " " << dd << " " << rr << G4endl;
    G4cout << "t1,Fctcos " << t1 << " " << Fctcos(t1, aa, bb, cc, dd, rr) << G4endl;
    G4cout << "t2,Fctcos " << t2 << " " << Fctcos(t2, aa, bb, cc, dd, rr) << G4endl;
  }
  G4double eps = 0.001;
  G4int ind1 = 10;
  G4double t = 0.0;
  G4int ier1;
  ier1 = Rtmi(&t, t1, t2, eps, ind1,
	      aa, bb, cc, dd, rr);
  if (verboseLevel > 1) {
    G4cout << "From Rtmi, ier1=" << ier1 << " t= " << t << G4endl;
    G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << G4endl;
  }
  if (ier1 != 0) t = 0.25*(3.*t1 + t2);
  if (verboseLevel > 1) {
      G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << 
              G4endl;
  }
  return t;
}

// The following is a "translation" of a root-finding routine
// from GEANT3.21/GHEISHA.  Some of the labelled block structure has
// been retained for clarity.  This routine will not be needed after
// the planned revisions to DoIt().

G4int
G4HadronElastic::Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
		      G4int iend, 
		      G4double aa, G4double bb, G4double cc, G4double dd, 
		      G4double rr)
{
   G4int ier = 0;
   G4double xl = xli;
   G4double xr = xri;
   *x = xl;
   G4double tol = *x;
   G4double f = Fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;
   G4double fl, fr;
   fl = f;
   *x = xr;
   tol = *x;
   f = Fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;
   fr = f;

// Error return in case of wrong input data
   if (fl*fr >= 0.) {
      ier = 2;
      return ier;
   }

// Basic assumption fl*fr less than 0 is satisfied.
// Generate tolerance for function values.
   G4int i = 0;
   G4double tolf = 100.*eps;

// Start iteration loop
label4:
   i++;

// Start bisection loop
   for (G4int k = 1; k <= iend; k++) {
      *x = 0.5*(xl + xr);
      tol = *x;
      f = Fctcos(tol, aa, bb, cc, dd, rr);
      if (f == 0.) return 0;
      if (f*fr < 0.) {      // Interchange xl and xr in order to get the
         tol = xl;          // same Sign in f and fr
         xl = xr;
         xr = tol;
         tol = fl;
         fl = fr;
         fr = tol;
      }
      tol = f - fl;
      G4double a = f*tol;
      a = a + a;
      if (a < fr*(fr - fl) && i <= iend) goto label17;
      xr = *x;
      fr = f;

// Test on satisfactory accuracy in bisection loop
      tol = eps;
      a = std::abs(xr);
      if (a > 1.) tol = tol*a;
      if (std::abs(xr - xl) <= tol && std::abs(fr - fl) <= tolf) goto label14;
   }
// End of bisection loop

// No convergence after iend iteration steps followed by iend
// successive steps of bisection or steadily increasing function
// values at right bounds.  Error return.
   ier = 1;

label14:
   if (std::abs(fr) > std::abs(fl)) {
      *x = xl;
      f = fl;
   }
   return ier;

// Computation of iterated x-value by inverse parabolic interp
label17:
   G4double a = fr - f;
   G4double dx = (*x - xl)*fl*(1. + f*(a - tol)/(a*(fr - fl)))/tol;
   G4double xm = *x;
   G4double fm = f;
   *x = xl - dx;
   tol = *x;
   f = Fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;

// Test on satisfactory accuracy in iteration loop
   tol = eps;
   a = std::abs(*x);
   if (a > 1) tol = tol*a;
   if (std::abs(dx) <= tol && std::abs(f) <= tolf) return ier;

// Preparation of next bisection loop
   if (f*fl < 0.) {
      xr = *x;
      fr = f;
   }
   else {
      xl = *x;
      fl = f;
      xr = xm;
      fr = fm;
   }
   goto label4;
}

// Test function for root-finder
G4double
G4HadronElastic::Fctcos(G4double t, 
			G4double aa, G4double bb, G4double cc, G4double dd, 
			G4double rr)
{
   const G4double expxl = -82.;
   const G4double expxu = 82.;

   G4double test1 = -bb*t;
   if (test1 > expxu) test1 = expxu;
   if (test1 < expxl) test1 = expxl;

   G4double test2 = -dd*t;
   if (test2 > expxu) test2 = expxu;
   if (test2 < expxl) test2 = expxl;

   return aa*std::exp(test1) + cc*std::exp(test2) - rr;
}


