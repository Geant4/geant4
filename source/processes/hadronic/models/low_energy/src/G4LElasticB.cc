//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// Physics model class G4LElasticB (derived from G4LElastic)
//
// G4 Model: Low-energy Elastic scattering with 4-momentum balance
// F.W. Jones, TRIUMF, 04-JUN-96
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Elastic class
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
//

#include "G4LElasticB.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4HadFinalState*
G4LElasticB::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
   if(getenv("debug_LElastic")) verboseLevel = 5;
   theParticleChange.Clear();
   const G4HadProjectile* aParticle = &aTrack;
   G4double atno2 = targetNucleus.GetN();
   G4double zTarget = targetNucleus.GetZ();
   theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
   theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());

// Elastic scattering off Hydrogen

   G4DynamicParticle* aSecondary = 0;
   if (atno2 < 1.5) {
      const G4ParticleDefinition* aParticleType = aParticle->GetDefinition();
      if (aParticleType == G4PionPlus::PionPlus())
         aSecondary = LightMedia.PionPlusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4PionMinus::PionMinus())
         aSecondary = LightMedia.PionMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4KaonPlus::KaonPlus())
         aSecondary = LightMedia.KaonPlusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4KaonZeroShort::KaonZeroShort())
         aSecondary = LightMedia.KaonZeroShortExchange(aParticle,targetNucleus);
      else if (aParticleType == G4KaonZeroLong::KaonZeroLong())
         aSecondary = LightMedia.KaonZeroLongExchange(aParticle, targetNucleus);
      else if (aParticleType == G4KaonMinus::KaonMinus())
         aSecondary = LightMedia.KaonMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4Proton::Proton())
         aSecondary = LightMedia.ProtonExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiProton::AntiProton())
         aSecondary = LightMedia.AntiProtonExchange(aParticle, targetNucleus);
      else if (aParticleType == G4Neutron::Neutron())
         aSecondary = LightMedia.NeutronExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiNeutron::AntiNeutron())
         aSecondary = LightMedia.AntiNeutronExchange(aParticle, targetNucleus);
      else if (aParticleType == G4Lambda::Lambda())
         aSecondary = LightMedia.LambdaExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiLambda::AntiLambda())
         aSecondary = LightMedia.AntiLambdaExchange(aParticle, targetNucleus);
      else if (aParticleType == G4SigmaPlus::SigmaPlus())
         aSecondary = LightMedia.SigmaPlusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4SigmaMinus::SigmaMinus())
         aSecondary = LightMedia.SigmaMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiSigmaPlus::AntiSigmaPlus())
         aSecondary = LightMedia.AntiSigmaPlusExchange(aParticle,targetNucleus);
      else if (aParticleType == G4AntiSigmaMinus::AntiSigmaMinus())
         aSecondary= LightMedia.AntiSigmaMinusExchange(aParticle,targetNucleus);
      else if (aParticleType == G4XiZero::XiZero())
         aSecondary = LightMedia.XiZeroExchange(aParticle, targetNucleus);
      else if (aParticleType == G4XiMinus::XiMinus())
         aSecondary = LightMedia.XiMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiXiZero::AntiXiZero())
         aSecondary = LightMedia.AntiXiZeroExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiXiMinus::AntiXiMinus())
         aSecondary = LightMedia.AntiXiMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4OmegaMinus::OmegaMinus())
         aSecondary = LightMedia.OmegaMinusExchange(aParticle, targetNucleus);
      else if (aParticleType == G4AntiOmegaMinus::AntiOmegaMinus())
         aSecondary= LightMedia.AntiOmegaMinusExchange(aParticle,targetNucleus);
      else if (aParticleType == G4KaonPlus::KaonPlus())
         aSecondary = LightMedia.KaonPlusExchange(aParticle, targetNucleus);
   }

   // Has a charge or strangeness exchange occurred?
   if (aSecondary) {
      aSecondary->SetMomentum(aParticle->Get4Momentum().vect());
      theParticleChange.SetStatusChange(stopAndKill);
      theParticleChange.AddSecondary(aSecondary);
   }
   // G4cout << "Entering elastic scattering 1"<<G4endl;

   G4double p = aParticle->GetTotalMomentum()/GeV;
   if (verboseLevel > 1)
      G4cout << "G4LElasticB::DoIt: Incident particle p=" << p << " GeV" << G4endl;

   if (p < 0.01) return &theParticleChange;

   // G4cout << "Entering elastic scattering 2"<<G4endl;
   // Compute the direction of elastic scattering.
   // It is planned to replace this code with a method based on
   // parameterized functions and a Monte Carlo method to invert the CDF.

   G4double ran = G4UniformRand();
   G4double aa, bb, cc, dd, rr;
   if (atno2 <= 62.) {
      aa = std::pow(atno2, 1.63);
      bb = 14.5*std::pow(atno2, 0.66);
      cc = 1.4*std::pow(atno2, 0.33);
      dd = 10.;
   }
   else {
      aa = std::pow(atno2, 1.33);
      bb = 60.*std::pow(atno2, 0.33);
      cc = 0.4*std::pow(atno2, 0.40);
      dd = 10.;
   }
   aa = aa/bb;
   cc = cc/dd;
   rr = (aa + cc)*ran;
   if (verboseLevel > 1) {
      G4cout << "DoIt: aa,bb,cc,dd,rr" << G4endl;
      G4cout << aa << " " << bb << " " << cc << " " << dd << " " << rr << G4endl;
   }
   G4double t1 = -std::log(ran)/bb;
   G4double t2 = -std::log(ran)/dd;
   if (verboseLevel > 1) {
      G4cout << "t1,Fctcos " << t1 << " " << Fctcos(t1, aa, bb, cc, dd, rr) << 
              G4endl;
      G4cout << "t2,Fctcos " << t2 << " " << Fctcos(t2, aa, bb, cc, dd, rr) << 
              G4endl;
   }
   G4double eps = 0.001;
   G4int ind1 = 10;
   G4double t = 0.0;
   G4int ier1;
   ier1 = Rtmi(&t, t1, t2, eps, ind1,
               aa, bb, cc, dd, rr);
   if (verboseLevel > 1) {
      G4cout << "From Rtmi, ier1=" << ier1 << G4endl;
      G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << 
              G4endl;
   }
   if (ier1 != 0) t = 0.25*(3.*t1 + t2);
   if (verboseLevel > 1) {
      G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << 
              G4endl;
   }

   // Scattered particle referred to axis of incident particle
   G4double m1=aParticle->GetDefinition()->GetPDGMass();
   if(aSecondary) m1 = aSecondary->GetMass();

   G4int Z=static_cast<G4int>(zTarget+.5);
   G4int A=static_cast<G4int>(atno2);
   G4ParticleDefinition * theDef = 0;
   if(Z == 1 && A == 1) {
     theDef = G4Proton::Proton();
   } else { 
     if(G4UniformRand()<atno2-A) A++;
     theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);
   }
   G4double m2 = theDef->GetPDGMass();
   G4LorentzVector lv1 = aParticle->Get4Momentum();
   G4LorentzVector lv0(0.0,0.0,0.0,m2);
   
   G4LorentzVector lv  = lv0 + lv1;
   G4ThreeVector bst = lv.boostVector();
   lv1.boost(-bst);
   lv0.boost(-bst);
   G4ThreeVector p1 = lv1.vect();
   G4double ptot = p1.mag();
   G4double ptotgev = ptot/GeV;
   if (ptotgev < 0.01) return &theParticleChange;

   // Sampling in CM system
   G4double phi  = G4UniformRand()*twopi;
   G4double cost = 1. - 0.5*t/(ptotgev*ptotgev);
   G4double sint = 0.0;
   if(cost > 1.0) cost = 1.0;
   else {
     if(cost < -1.0) cost = -1.0 + 2.0*G4UniformRand();
     sint = std::sqrt((1.0-cost)*(1.0+cost));
   }
   // G4cout << "Entering elastic scattering 3"<<G4endl;
   if (verboseLevel > 1) 
     G4cout << "cos(t)=" << cost << " sin(t)=" << sint << G4endl;

   G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
   p1 = p1.unit();
   v1.rotateUz(p1);
   v1 *= ptot;
   G4LorentzVector nlv11(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));
   G4LorentzVector nlv01 = lv0 + lv1 - nlv11;
   nlv01.boost(bst);
   nlv11.boost(bst); 
   //G4cout << " P0= "<< nlv01 << "   P1= "<< nlv11<<" m= " << m1 <<G4endl;

   if (aSecondary) {
     aSecondary->SetMomentum(nlv11.vect());
   } else {
     theParticleChange.SetMomentumChange(nlv11.vect().unit());
     theParticleChange.SetEnergyChange(std::max(nlv11.e() - m1,0.0));
   }
   G4DynamicParticle * aSec = new G4DynamicParticle(theDef, nlv01);
   theParticleChange.AddSecondary(aSec);

   //   G4cout << " ion info "<<atno2 << " "<<A<<" "<<Z<<" "<<zTarget<<G4endl;
   return &theParticleChange;
}

// The following is a "translation" of a root-finding routine
// from GEANT3.21/GHEISHA.  Some of the labelled block structure has
// been retained for clarity.  This routine will not be needed after
// the planned revisions to DoIt().

G4int
G4LElasticB::Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
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
G4LElasticB::Fctcos(G4double t, 
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


void
G4LElasticB::Defs1(G4double p, G4double px, G4double py, G4double pz, 
                  G4double pxinc, G4double pyinc, G4double pzinc, 
                  G4double* pxnew, G4double* pynew, G4double* pznew)
{
// Transform scattered particle to reflect direction of incident particle
   G4double pt2 = pxinc*pxinc + pyinc*pyinc;
   if (pt2 > 0.) {
      G4double cost = pzinc/p;
      G4double sint1 = std::sqrt(std::abs((1. - cost )*(1.+cost)));
      G4double sint2 = std::sqrt(pt2)/p;
      G4double sint = 0.5*(sint1 + sint2);
      G4double ph = pi*0.5;
      if (pyinc < 0.) ph = pi*1.5;
      if (std::abs(pxinc) > 1.e-6) ph = std::atan2(pyinc, pxinc);
      G4double cosp = std::cos(ph);
      G4double sinp = std::sin(ph);
      if (verboseLevel > 1) {
         G4cout << "cost sint " << cost << " " << sint << G4endl;
         G4cout << "cosp sinp " << cosp << " " << sinp << G4endl;
      }
      *pxnew = cost*cosp*px - sinp*py + sint*cosp*pz;
      *pynew = cost*sinp*px + cosp*py + sint*sinp*pz;
      *pznew =     -sint*px                 +cost*pz;
   }
   else {
       G4double cost=pzinc/p;
       *pxnew = cost*px;
       *pynew = py;
       *pznew = cost*pz;
   }
}
