// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPMadlandNixSpectrum.hh,v 1.4 1999-12-15 14:53:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPMadlandNixSpectrum_h
#define G4NeutronHPMadlandNixSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include <math.h>
#include "G4VNeutronHPEDis.hh"

//     #include <nag.h> @
//     #include <nags.h> @


// we will need a List of these .... one per term.

class G4NeutronHPMadlandNixSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPMadlandNixSpectrum()
  {
    expm1 = exp(-1.);
  }
  ~G4NeutronHPMadlandNixSpectrum()
  {
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile);
    aDataFile>> theAvarageKineticPerNucleonForLightFragments;
    theAvarageKineticPerNucleonForLightFragments*=eV;
    aDataFile>> theAvarageKineticPerNucleonForHeavyFragments;
    theAvarageKineticPerNucleonForHeavyFragments*=eV;
    theMaxTemp.Init(aDataFile);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy);
    
  private:
  
  G4double Madland(G4double aSecEnergy, G4double tm);
  
  inline G4double FissionIntegral(G4double tm, G4double anEnergy)
  {
    return 0.5*(  GIntegral(tm, anEnergy, theAvarageKineticPerNucleonForLightFragments) 
                       +GIntegral(tm, anEnergy, theAvarageKineticPerNucleonForHeavyFragments) );
  }
  
  G4double GIntegral(G4double tm, G4double anEnergy, G4double aMean);
  
  inline G4double Gamma05(G4double aValue)
  {
    G4double result;
    // gamma(1.2,x*X) = sqrt(pi)*Erf(x)
    G4double x = sqrt(aValue);
    G4double t = 1./(1+0.47047*x);
    result = 1- (0.3480242*t - 0.0958798*t*t + 0.7478556*t*t*t)*exp(-aValue); // @ check
    result *= sqrt(pi);
    return result;
  }
  
  inline G4double Gamma15(G4double aValue)
  {
    G4double result;
    // gamma(a+1, x) = a*gamma(a,x)-x**a*exp(-x)
    result = 0.5*Gamma05(aValue) - sqrt(aValue)*exp(-aValue); // @ check
    return result;
  }
  
  inline G4double Gamma25(G4double aValue)
  {
    G4double result;
    result = 1.5*Gamma15(aValue) - pow(aValue,1.5)*exp(aValue); // @ check
    return result;
  }
  
  inline G4double E1(G4double aValue)
  {
  // good only for rather low aValue @@@ replace by the corresponding NAG function for the
  // exponential integral. (<5 seems ok.
    G4double gamma = 0.577216;
    G4double precision = 0.000001;
    G4double result =-gamma - log(aValue);
    G4double term = -aValue;
    G4double last;
    G4int count = 1;
    result -= term;
    for(;;)
    {
      count++;
      last = result;
      term = -term*aValue*(count-1)/(count*count);
      result -=term;
      if(fabs(term)/fabs(result)<precision) break;
    }
//    NagError *fail; @
//    result = nag_exp_integral(aValue, fail); @
    return result;
  }
  
  private:
  
  G4double expm1;
  
  private: 
  
  G4NeutronHPVector theFractionalProb;
  
  G4double theAvarageKineticPerNucleonForLightFragments;
  G4double theAvarageKineticPerNucleonForHeavyFragments;

  G4NeutronHPVector theMaxTemp;
  
};

#endif
