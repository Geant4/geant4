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
// $Id: G4eIonisationSpectrum.cc,v 1.9 2001-11-29 19:01:37 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eIonisationSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 29 September 2001
//
// Modifications: 
// 10.10.2001 MGP           Revision to improve code quality and 
//                          consistency with design
// 02.11.2001 VI            Optimize sampling of energy 
// 29.11.2001 VI            New parametrisation 
//
// -------------------------------------------------------------------
//

#include "G4eIonisationSpectrum.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4eIonisationParameters.hh"
#include "G4DataVector.hh"
#include "Randomize.hh"


G4eIonisationSpectrum::G4eIonisationSpectrum():G4VEnergySpectrum(),
  lowestE(0.1*eV),
  verbose(0)
{
  theParam = new G4eIonisationParameters();
}


G4eIonisationSpectrum::~G4eIonisationSpectrum() 
{
  delete theParam;
  theParam = 0;
}


G4double G4eIonisationSpectrum::Probability(G4int Z, 
  	    				    G4double tMin, 
					    G4double tMax, 
					    G4double e,
					    G4int shell,
				      const G4ParticleDefinition* part) const
{
  // Please comment what Probability does and what are the three 
  // functions mentioned below
  // Describe the algorithms used

  G4double eMax = MaxEnergyOfSecondaries(e);
  G4double t0 = G4std::max(tMin, lowestE);
  G4double tm = G4std::min(tMax, eMax);
  if(t0 >= tm) return 0.0;

  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  G4double x1 = (t0 + bindingEnergy)/(e + bindingEnergy);
  G4double x2 = G4std::min(0.5,(tm + bindingEnergy)/(e + bindingEnergy));

  if(verbose > 1) {
    G4cout << "G4eIonisationSpectrum::Probability: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; x1= " << x1 
           << "; x2= " << x2 
           << G4endl;
  }

  G4int iMax = 7;
  G4DataVector p;
  p.clear();

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }

  G4double g = (e + bindingEnergy)/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));
    
  G4double val = IntSpectrum(x1, x2, p);
  G4double x0  = (lowestE + bindingEnergy)/(e + bindingEnergy);
  G4double nor = IntSpectrum(x0, 0.5, p);
  
  if(verbose > 1) {
    G4cout << "tcut= " << tMin 
           << "; tMax= " << tMax 
           << "; x0= " << x0 
           << "; x1= " << x1 
           << "; x2= " << x2 
           << "; val= " << val 
           << "; nor= " << nor 
           << "; sum= " << p[0] 
           << "; a= " << p[1] 
           << "; b= " << p[2] 
           << "; c= " << p[3] 
           << G4endl;
  }

  if(nor > 0.0) val /= nor;
  else          val  = 0.0;
  if(val < 0.0) val  = 0.0;

  return val; 
}


G4double G4eIonisationSpectrum::AverageEnergy(G4int Z,
	        			      G4double tMin, 
					      G4double tMax, 
					      G4double e,
					      G4int shell,
				        const G4ParticleDefinition* part) const
{
  // Please comment what AverageEnergy does and what are the three 
  // functions mentioned below
  // Describe the algorithms used

  G4double eMax = MaxEnergyOfSecondaries(e);
  G4double t0 = G4std::max(tMin, lowestE);
  G4double tm = G4std::min(tMax, eMax);
  if(t0 >= tm) return 0.0;

  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  G4double x1 = (t0 + bindingEnergy)/(e + bindingEnergy);
  G4double x2 = G4std::min(0.5,(tm + bindingEnergy)/(e + bindingEnergy));

  if(verbose > 1) {
    G4cout << "G4eIonisationSpectrum::AverageEnergy: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; bindingE(keV)= " << bindingEnergy/keV
           << "; x1= " << x1 
           << "; x2= " << x2 
           << G4endl;
  }

  G4int iMax = 7;
  G4DataVector p;
  p.clear();

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }

  G4double g = (e + bindingEnergy)/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));
    
  G4double val = AverageValue(x1, x2, p);
  G4double x0  = (lowestE + bindingEnergy)/(e + bindingEnergy);
  G4double nor = IntSpectrum(x0, 0.5, p);
  val *= (e + bindingEnergy);

  if(verbose > 1) {
    G4cout << "tcut(MeV)= " << tMin/MeV 
           << "; tMax(MeV)= " << tMax/MeV 
           << "; x0= " << x0 
           << "; x1= " << x1 
           << "; x2= " << x2 
           << "; val= " << val 
           << "; nor= " << nor 
           << "; sum= " << p[0] 
           << "; a= " << p[1] 
           << "; b= " << p[2] 
           << "; c= " << p[3] 
           << G4endl;
  }

  if(nor > 0.0) val /= nor;
  else          val  = 0.0;
  if(val < 0.0) val  = 0.0;

  return val; 
}


G4double G4eIonisationSpectrum::SampleEnergy(G4int Z,
					     G4double tMin, 
					     G4double tMax, 
				             G4double e,
				             G4int shell,
				       const G4ParticleDefinition* part) const
{
  // Please comment what SampleEnergy does
  G4double tDelta = 0.0;
  G4double t0 = G4std::max(tMin, lowestE);
  G4double tm = G4std::min(tMax, MaxEnergyOfSecondaries(e));
  if(t0 > tm) return tDelta;

  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  G4double x1 = (t0 + bindingEnergy)/(e + bindingEnergy);
  G4double x2 = G4std::min(0.5,(tm + bindingEnergy)/(e + bindingEnergy));

  if(verbose > 1) {
    G4cout << "G4eIonisationSpectrum::SampleEnergy: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << G4endl;
  }

  // Access parameters
  G4int iMax = 7;
  G4DataVector p;
  p.clear();

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }

  G4double g = (e + bindingEnergy)/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));


  G4double aria1 = 0.0;
  G4double a1 = G4std::min(x1,p[6]);
  G4double a2 = G4std::min(x2,p[6]);
  if(a1 < a2) aria1 = IntSpectrum(a1, a2, p);
  G4double aria2 = 0.0;
  G4double a3 = G4std::max(x1,p[6]);
  G4double a4 = G4std::max(x2,p[6]);
  if(a3 < a4) aria2 = IntSpectrum(a3, a4, p);

  G4double aria = (aria1 + aria2)*G4UniformRand();
  G4double amaj, fun, q, x;

  //======= First aria to sample =====

  if(aria <= aria1) { 

    amaj = p[4];
    a1 = 1./a1;
    a2 = 1./a2;

  //======= Second aria to sample =====

  } else {

    amaj = p[5];
    a1 = 1./a3;
    a2 = 1./a4;
  }
  amaj *= 1.1;

  do {

    x = 1./(a2 + G4UniformRand()*(a1 - a2));
    fun  = Function(x, p);

    if(fun > amaj) {
          G4cout << "WARNING in G4eIonisationSpectrum::SampleEnergy:" 
                 << " Majoranta " << amaj 
                 << " < " << fun
                 << G4endl;
    }

    q = amaj*G4UniformRand();

  } while (q >= fun);

  p.clear();

  tDelta = x*(e + bindingEnergy) - bindingEnergy;

  return tDelta; 
}


G4double G4eIonisationSpectrum::IntSpectrum(G4double xMin, 
					    G4double xMax,
				      const G4DataVector& p) const
{
  // Please comment what IntSpectrum does
  G4double x1 = 1./xMin;
  G4double x2 = 1./xMax;
  G4double x  = x1 - x2 - p[7]*log(xMax/xMin) + (1. - p[7])*(xMax - xMin)
              + 1./(1. - xMax) - 1./(1. - xMin) 
              + p[7]*log((1. - xMax)/(1. - xMin))
              + 0.5*p[1]*p[3]*(x1*x1 - x2*x2);

  if(x < 0.0) x = 0.0;
  return x;
} 


G4double G4eIonisationSpectrum::AverageValue(G4double xMin, 
				             G4double xMax,
				       const G4DataVector& p) const
{

  G4double x1 = 1.;
  G4double x2 = 1.;
  G4double x  = log(xMax/xMin) 
              + 0.5*(1. - p[7])*(xMax*xMax - xMin*xMin)
              + 1./(1. - xMax) - 1./(1. - xMin) 
              + (1. + p[7])*log((1. - xMax)/(1. - xMin))
              + p[1]*p[3]*(1./xMin - 1./xMax);

  if(x < 0.0) x = 0.0;
  return x;
} 


G4double G4eIonisationSpectrum::Function(G4double x, 
				   const G4DataVector& p) const
{
  // Please comment what Function does


  G4double x1 = 1.0;
  G4double f  = 1.0 - p[7]*x + x*x*(1.0 - p[7]
              + (1.0/(1.0 - x) - p[7])/(1.0 - x) )
              + p[1]*p[3]/x;

  if(f < 0.0) f = 0.0;
  return f;
} 

G4double G4eIonisationSpectrum::Excitation(G4int Z, G4double e) const 
{
  return theParam->Excitation(Z, e);
}

void G4eIonisationSpectrum::PrintData() const 
{
  theParam->PrintData();
}





