#ifndef G4Rutherford_hh
#define G4Rutherford_hh

// J.P. Wellisch, X-mas 2002.
#include "G4HadronicInteraction.hh"
#include "G4ParticleChange.hh"

class G4Rutherford : public G4HadronicInteraction
{
   public:
   G4Rutherford(G4double cut=0.9999) 
   {
     game = false;
     theMaxCosTh=cut;
   }
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);
				    
   G4double Apply(const G4ParticleDefinition * aP, G4Nucleus& targetNucleus);
   private:

   double Ruther(double r, double cth)
   {
  //
  // dsig/dOm = (z1*z2*e^2/(16*pi*eps0*Elab))^2
  //           * 4./sin(th)^4
  //           *(sqrt( 1-(m1/m2)^2*sin(th)^2 ) + cos(th))^2
  //           / sqrt( 1-(m1/m2)^2*sin(th)^2 )
  //
     double result = 0;
     double cth2 = cth*cth;
     result = 1./(1-cth2)/(1-cth2);
     double a=1-r*r;
     double fac = sqrt(a+r*r*cth2);
     result *= (fac+cth)*(fac+cth);
     result /= fac;
     return result;
   }

   G4ParticleChange theResult;
   G4bool game;
   G4double theMaxCosTh;
};

#endif

