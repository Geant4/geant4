#ifndef G4Rutherford_hh
#define G4Rutherford_hh

// J.P. Wellisch, X-mas 2002.
#include "G4HadronicInteraction.hh"
#include "G4ParticleChange.hh"
#include "G4VCrossSectionDataSet.hh"

class G4Rutherford : public G4HadronicInteraction, public G4VCrossSectionDataSet
{
   public:
   G4Rutherford(G4double cut=0.99) 
   {
     game = false;
     theMaxCosTh=cut;
   }
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);
				    
   G4double Apply(const G4ParticleDefinition * aP, G4Nucleus& targetNucleus);
   virtual G4bool IsApplicable(const G4DynamicParticle* aP, const G4Element*aE)
   {
     G4bool result=false;
     if(aP->GetDefinition()->GetBaryonNumber()>1.5)
     {
       result = true;
     }
     return result;
   }

   virtual G4double GetCrossSection(const G4DynamicParticle* aP,
                                   const G4Element* aE,
                                   G4double aTemperature = 0.)
   {
     G4Nucleus aNuc(aE->GetN(), aE->GetZ());
     return XSec(aP->GetKineticEnergy(), aP->GetDefinition(), aNuc);
   }
   
   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&aP){}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&aP){};
   G4double XSec(G4double eLab, const G4ParticleDefinition * aP, G4Nucleus& targetNucleus)
   {
     prepare(aP, targetNucleus);
     G4double z1=targetNucleus.GetZ();
     G4double z2=aP->GetPDGCharge();
     G4double result = z1*z2*electron_charge*electron_charge;
     result /= 16.*pi*epsilon0*eLab;
     result *= result;
     result *= 2.*pi;
     result *= -total;
     G4cout << "total = "<<total<<G4endl;
     return result;
   }
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
   
   void prepare(const G4ParticleDefinition* aP, G4Nucleus& targetNucleus);

   G4ParticleChange theResult;
   G4bool game;
   G4double theMaxCosTh;
   
   vector<G4double> theRuther;
   vector<G4double> theValue;
   vector<G4double> integral;
   G4double total;

   
};

#endif

