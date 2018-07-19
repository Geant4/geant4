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
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPMadlandNixSpectrum_h
#define G4ParticleHPMadlandNixSpectrum_h 1

#include <fstream>
#include <cmath>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4ParticleHPVector.hh"
#include "G4VParticleHPEDis.hh"

//     #include <nag.h> @
//     #include <nags.h> @


// we will need a List of these .... one per term.

class G4ParticleHPMadlandNixSpectrum : public G4VParticleHPEDis
{
  public:
  G4ParticleHPMadlandNixSpectrum()
  {
    expm1 = G4Exp(-1.);
    theAvarageKineticPerNucleonForLightFragments = 0.0;
    theAvarageKineticPerNucleonForHeavyFragments = 0.0;
  }
  ~G4ParticleHPMadlandNixSpectrum()
  {
  }
  
  inline void Init(std::istream & aDataFile)
  {
    theFractionalProb.Init(aDataFile);
    aDataFile>> theAvarageKineticPerNucleonForLightFragments;
    theAvarageKineticPerNucleonForLightFragments*=CLHEP::eV;
    aDataFile>> theAvarageKineticPerNucleonForHeavyFragments;
    theAvarageKineticPerNucleonForHeavyFragments*=CLHEP::eV;
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
    // gamma(1.2,x*X) = std::sqrt(CLHEP::pi)*Erf(x)
    G4double x = std::sqrt(aValue);
    G4double t = 1./(1+0.47047*x);
    result = 1- (0.3480242*t - 0.0958798*t*t + 0.7478556*t*t*t)*G4Exp(-aValue); // @ check
    result *= std::sqrt(CLHEP::pi);
    return result;
  }
  
  inline G4double Gamma15(G4double aValue)
  {
    G4double result;
    // gamma(a+1, x) = a*gamma(a,x)-x**a*std::exp(-x)
    result = 0.5*Gamma05(aValue) - std::sqrt(aValue)*G4Exp(-aValue); // @ check
    return result;
  }
  
  inline G4double Gamma25(G4double aValue)
  {
    G4double result;
    result = 1.5*Gamma15(aValue) - G4Pow::GetInstance()->powA(aValue,1.5)*G4Exp(aValue); // @ check
    return result;
  }
  
  inline G4double E1(G4double aValue)
  {
  // good only for rather low aValue @@@ replace by the corresponding NAG function for the
  // exponential integral. (<5 seems ok.
    G4double gamma = 0.577216;
    G4double precision = 0.000001;
    G4double result =-gamma - G4Log(aValue);
    G4double term = -aValue;
    //110527TKDB  Unnessary codes, Detected by gcc4.6 compiler 
    //G4double last;
    G4int count = 1;
    result -= term;
    for(;;)
    {
      count++;
      //110527TKDB  Unnessary codes, Detected by gcc4.6 compiler 
      //last = result;
      term = -term*aValue*(count-1)/(count*count);
      result -=term;
      if(std::fabs(term)/std::fabs(result)<precision) break;
    }
//    NagError *fail; @
//    result = nag_exp_integral(aValue, fail); @
    return result;
  }
  
  private:
  
  G4double expm1;
  
  private: 
  
  G4ParticleHPVector theFractionalProb;
  
  G4double theAvarageKineticPerNucleonForLightFragments;
  G4double theAvarageKineticPerNucleonForHeavyFragments;

  G4ParticleHPVector theMaxTemp;
  
};

#endif
