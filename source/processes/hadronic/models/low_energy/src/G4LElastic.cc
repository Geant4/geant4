// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LElastic.cc,v 1.2 1999-05-25 00:36:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Physics model class G4LElastic
//
//
// G4 Model: Low-energy Elastic scattering
// F.W. Jones, TRIUMF, 04-JUN-96
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Elastic class
//
// 25-JUN-98 FWJ: replaced missing Initialize for ParticleChange.
//

#include "globals.hh"
#include "G4LElastic.hh"
#include "Randomize.hh"


G4VParticleChange*
G4LElastic::ApplyYourself(const G4Track& aTrack, G4Nucleus& targetNucleus)
{
   theParticleChange.Initialize(aTrack);

   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4double Z = targetNucleus.GetZ();
   G4double atno2 = targetNucleus.GetN();

// Elastic scattering off Hydrogen

   G4DynamicParticle* aSecondary = 0;
   if (atno2 < 1.5) {
      G4ParticleDefinition* aParticleType = aParticle->GetDefinition();
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
      aSecondary->SetMomentum(aParticle->GetMomentum());
      theParticleChange.SetStatusChange(fStopAndKill);
      theParticleChange.AddSecondary(aSecondary);
   }

   G4double p = aParticle->GetTotalMomentum()/GeV;
   if (verboseLevel > 1)
      G4cout << "G4LElastic::DoIt: Incident particle p=" << p << " GeV" << endl;

   if (p < 0.01) return &theParticleChange;

// Compute the direction of elastic scattering.
// It is planned to replace this code with a method based on
// parameterized functions and a Monte Carlo method to invert the CDF.

   G4double ran = G4UniformRand();
   G4double aa, bb, cc, dd, rr;
   if (atno2 <= 62.) {
      aa = pow(atno2, 1.63);
      bb = 14.5*pow(atno2, 0.66);
      cc = 1.4*pow(atno2, 0.33);
      dd = 10.;
   }
   else {
      aa = pow(atno2, 1.33);
      bb = 60.*pow(atno2, 0.33);
      cc = 0.4*pow(atno2, 0.40);
      dd = 10.;
   }
   aa = aa/bb;
   cc = cc/dd;
   rr = (aa + cc)*ran;
   if (verboseLevel > 1) {
      G4cout << "DoIt: aa,bb,cc,dd,rr" << endl;
      G4cout << aa << " " << bb << " " << cc << " " << dd << " " << rr << endl;
   }
   G4double t1 = -log(ran)/bb;
   G4double t2 = -log(ran)/dd;
   if (verboseLevel > 1) {
      G4cout << "log(FLT_MAX)=" << log(FLT_MAX) << endl;
      G4cout << "t1,Fctcos " << t1 << " " << Fctcos(t1, aa, bb, cc, dd, rr) << 
              endl;
      G4cout << "t2,Fctcos " << t2 << " " << Fctcos(t2, aa, bb, cc, dd, rr) << 
              endl;
   }
   G4double eps = 0.001;
   G4int ind1 = 10;
   G4double t, val;
   G4int ier1;
   ier1 = Rtmi(&t, t1, t2, eps, ind1,
               aa, bb, cc, dd, rr);
   if (verboseLevel > 1) {
      G4cout << "From Rtmi, ier1=" << ier1 << endl;
      G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << 
              endl;
   }
   if (ier1 != 0) t = 0.25*(3.*t1 + t2);
   if (verboseLevel > 1) {
      G4cout << "t, Fctcos " << t << " " << Fctcos(t, aa, bb, cc, dd, rr) << 
              endl;
   }
   G4double phi = G4UniformRand()*twopi;
   rr = 0.5*t/(p*p);
   if (rr > 1.) rr = 0.;
   if (verboseLevel > 1)
      G4cout << "rr=" << rr << endl;
   G4double cost = 1. - rr;
   G4double sint = sqrt(max(rr*(2. - rr), 0.));
   if (sint == 0.) return &theParticleChange;
   if (verboseLevel > 1)
      G4cout << "cos(t)=" << cost << "  sin(t)=" << sint << endl;
// Scattered particle referred to axis of incident particle
   G4double px = p*sint*sin(phi);
   G4double py = p*sint*cos(phi);
   G4double pz = p*cost;
// Incident particle
   G4double pxinc = p*(aParticle->GetMomentumDirection().x());
   G4double pyinc = p*(aParticle->GetMomentumDirection().y());
   G4double pzinc = p*(aParticle->GetMomentumDirection().z());
   if (verboseLevel > 1) {
      G4cout << "NOM SCAT " << px << " " << py << " " << pz << endl;
      G4cout << "INCIDENT " << pxinc << " " << pyinc << " " << pzinc << endl;
   }

// Transform scattered particle to reflect direction of incident particle
   G4double pxnew, pynew, pznew;
   Defs1(p, px, py, pz, pxinc, pyinc, pzinc, &pxnew, &pynew, &pznew);
// Normalize:
   pxnew = pxnew/p;
   pynew = pynew/p;
   pznew = pznew/p;
   if (verboseLevel > 1) {
      G4cout << "DoIt: returning new momentum vector" << endl;
      G4cout << pxnew << " " << pynew << " " << pznew << endl;
   }

   if (aSecondary)
      aSecondary->SetMomentumDirection(pxnew, pynew, pznew);
   else
      theParticleChange.SetMomentumChange(pxnew, pynew, pznew);

   return &theParticleChange;
}


// The following is a "translation" of a root-finding routine
// from GEANT3.21/GHEISHA.  Some of the labelled block structure has
// been retained for clarity.  This routine will not be needed after
// the planned revisions to DoIt().

G4int
G4LElastic::Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
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
      a = abs(xr);
      if (a > 1.) tol = tol*a;
      if (abs(xr - xl) <= tol && abs(fr - fl) <= tolf) goto label14;
   }
// End of bisection loop

// No convergence after iend iteration steps followed by iend
// successive steps of bisection or steadily increasing function
// values at right bounds.  Error return.
   ier = 1;

label14:
   if (abs(fr) > abs(fl)) {
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
   a = abs(*x);
   if (a > 1) tol = tol*a;
   if (abs(dx) <= tol && abs(f) <= tolf) return ier;

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
G4LElastic::Fctcos(G4double t, 
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

   return aa*exp(test1) + cc*exp(test2) - rr;
}


void
G4LElastic::Defs1(G4double p, G4double px, G4double py, G4double pz, 
                  G4double pxinc, G4double pyinc, G4double pzinc, 
                  G4double* pxnew, G4double* pynew, G4double* pznew)
{
// Transform scattered particle to reflect direction of incident particle
   G4double pt2 = pxinc*pxinc + pyinc*pyinc;
   if (pt2 > 0.) {
      G4double cost = pzinc/p;
      G4double sint1 = sqrt(abs((1. - cost )*(1.+cost)));
      G4double sint2 = sqrt(pt2)/p;
      G4double sint = 0.5*(sint1 + sint2);
      G4double ph = pi*0.5;
      if (pyinc < 0.) ph = pi*1.5;
      if (abs(pxinc) > 1.e-6) ph = atan2(pyinc, pxinc);
      G4double cosp = cos(ph);
      G4double sinp = sin(ph);
      if (verboseLevel > 1) {
         G4cout << "cost sint " << cost << " " << sint << endl;
         G4cout << "cosp sinp " << cosp << " " << sinp << endl;
      }
      *pxnew = cost*cosp*px - sinp*py + sint*cosp*pz;
      *pynew = cost*sinp*px + cosp*py + sint*sinp*pz;
      *pznew =     -sint*px                 +cost*pz;
   }
   else {
      *pxnew = px;
      *pynew = py;
      *pznew = pz;
   }
}
