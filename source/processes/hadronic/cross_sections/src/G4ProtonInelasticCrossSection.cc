// By JPW, working, but to be cleaned up. @@@

#include "G4ProtonInelasticCrossSection.hh"
#include "globals.hh"

   G4double G4ProtonInelasticCrossSection::
   GetCrossSection(const G4DynamicParticle* aPart, const G4Element* anEle)
   {
      G4double atomicNumber = anEle->GetN();
      G4double nOfProtons = anEle->GetZ();
      return GetCrossSection(aPart->GetKineticEnergy(), atomicNumber, nOfProtons);
   }
   
   G4double G4ProtonInelasticCrossSection::
   GetCrossSection(G4double kineticEnergy, G4double atomicNumber, G4double nOfProtons)
   {
      G4double nOfNeutrons = atomicNumber-nOfProtons;
      kineticEnergy /=GeV;
      G4double a = atomicNumber;
      const G4double nuleonRadius=1.36E-15;
      const G4double pi=3.14159265;
      G4double fac=pi*nuleonRadius*nuleonRadius;
      G4double b0=2.247-0.915*(1-pow(a,-0.3333));
      G4double fac1=b0*(1-pow(a,-0.3333));
      G4double fac2=1.;
      if(nOfNeutrons>1.5) fac2=log((nOfNeutrons));
      G4double crossSection = 1E31*fac*fac2*(1+pow(a,0.3333)-fac1);
// high energy correction
      crossSection = (1-0.15*exp(-kineticEnergy))*crossSection/(1.00-0.0007*a);
// first try on low energies: rise
      G4double ff1= 0.70-0.002*a;  // slope of the drop at medium energies.
      G4double ff2= 1.00+1/a;  // start of the slope.
      G4double ff3= 0.8+18/a-0.002*a; // stephight
      fac=1-(1/(1+exp(-8*ff1*(log10(kineticEnergy)+1.37*ff2))));
      crossSection = crossSection*(1+ff3*fac);
// low energy return to zero
      ff1=1.-1/a-0.001*a; // slope of the rise
      ff2=1.17-2.7/a-0.0014*a; // start of the rise
      fac=-8.*ff1*(log10(kineticEnergy)+2.0*ff2);
      fac=1/(1+exp(fac));
      crossSection = crossSection*fac;
      return crossSection*millibarn;
    }
