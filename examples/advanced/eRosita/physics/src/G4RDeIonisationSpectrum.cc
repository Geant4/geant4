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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4RDeIonisationSpectrum
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
// 19.04.2002 VI            Add protection in case of energy below binding  
// 30.05.2002 VI            Update to 24-parameters data
// 11.07.2002 VI            Fix in integration over spectrum
//
// -------------------------------------------------------------------
//

#include "G4RDeIonisationSpectrum.hh"
#include "G4RDAtomicTransitionManager.hh"
#include "G4RDAtomicShell.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "Randomize.hh"


G4RDeIonisationSpectrum::G4RDeIonisationSpectrum():G4RDVEnergySpectrum(),
  lowestE(0.1*eV),
  factor(1.3),
  iMax(24),
  verbose(0)
{
  theParam = new G4RDeIonisationParameters();
}


G4RDeIonisationSpectrum::~G4RDeIonisationSpectrum() 
{
  delete theParam;
}


G4double G4RDeIonisationSpectrum::Probability(G4int Z, 
  	    				    G4double tMin, 
					    G4double tMax, 
					    G4double e,
					    G4int shell,
					    const G4ParticleDefinition* ) const
{
  // Please comment what Probability does and what are the three 
  // functions mentioned below
  // Describe the algorithms used

  G4double eMax = MaxEnergyOfSecondaries(e);
  G4double t0 = std::max(tMin, lowestE);
  G4double tm = std::min(tMax, eMax);
  if(t0 >= tm) return 0.0;

  G4double bindingEnergy = (G4RDAtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  if(e <= bindingEnergy) return 0.0;

  G4double energy = e + bindingEnergy;

  G4double x1 = std::min(0.5,(t0 + bindingEnergy)/energy);
  G4double x2 = std::min(0.5,(tm + bindingEnergy)/energy);

  if(verbose > 1 || (Z==4 && e>= 1.0 && e<= 0.0)) {
    G4cout << "G4RDeIonisationSpectrum::Probability: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; Eb(keV)= " << bindingEnergy/keV
           << "; x1= " << x1 
           << "; x2= " << x2 
           << G4endl;
     
  }

  G4DataVector p;

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    G4double x = theParam->Parameter(Z, shell, i, e);
    if(i<4) x /= energy; 
    p.push_back(x); 
  }

  if(p[3] > 0.5) p[3] = 0.5;
  
  G4double g = energy/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));
  
  p[iMax-1] = Function(p[3], p);

  if(e >= 1. && e <= 0. && Z == 4) p.push_back(0.0);

  
  G4double val = IntSpectrum(x1, x2, p);
  G4double x0  = (lowestE + bindingEnergy)/energy;
  G4double nor = IntSpectrum(x0, 0.5, p);
  
  if(verbose > 1 || (Z==4 && e>= 1.0 && e<= 0.0)) {
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
    if(shell == 1) G4cout << "============" << G4endl; 
  }

  p.clear();

  if(nor > 0.0) val /= nor;
  else          val  = 0.0;

  return val; 
}


G4double G4RDeIonisationSpectrum::AverageEnergy(G4int Z,
	        			      G4double tMin, 
					      G4double tMax, 
					      G4double e,
					      G4int shell,
					      const G4ParticleDefinition* ) const
{
  // Please comment what AverageEnergy does and what are the three 
  // functions mentioned below
  // Describe the algorithms used

  G4double eMax = MaxEnergyOfSecondaries(e);
  G4double t0 = std::max(tMin, lowestE);
  G4double tm = std::min(tMax, eMax);
  if(t0 >= tm) return 0.0;

  G4double bindingEnergy = (G4RDAtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  if(e <= bindingEnergy) return 0.0;

  G4double energy = e + bindingEnergy;

  G4double x1 = std::min(0.5,(t0 + bindingEnergy)/energy);
  G4double x2 = std::min(0.5,(tm + bindingEnergy)/energy);

  if(verbose > 1) {
    G4cout << "G4RDeIonisationSpectrum::AverageEnergy: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; bindingE(keV)= " << bindingEnergy/keV
           << "; x1= " << x1 
           << "; x2= " << x2 
           << G4endl;
  }

  G4DataVector p;

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    G4double x = theParam->Parameter(Z, shell, i, e);
    if(i<4) x /= energy; 
    p.push_back(x);
  }

  if(p[3] > 0.5) p[3] = 0.5;

  G4double g = energy/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));

  p[iMax-1] = Function(p[3], p);
    
  G4double val = AverageValue(x1, x2, p);
  G4double x0  = (lowestE + bindingEnergy)/energy;
  G4double nor = IntSpectrum(x0, 0.5, p);
  val *= energy;

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

  p.clear();

  if(nor > 0.0) val /= nor;
  else          val  = 0.0;

  return val; 
}


G4double G4RDeIonisationSpectrum::SampleEnergy(G4int Z,
					     G4double tMin, 
					     G4double tMax, 
				             G4double e,
				             G4int shell,
					     const G4ParticleDefinition* ) const
{
  // Please comment what SampleEnergy does
  G4double tDelta = 0.0;
  G4double t0 = std::max(tMin, lowestE);
  G4double tm = std::min(tMax, MaxEnergyOfSecondaries(e));
  if(t0 > tm) return tDelta;

  G4double bindingEnergy = (G4RDAtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();

  if(e <= bindingEnergy) return 0.0;

  G4double energy = e + bindingEnergy;

  G4double x1 = std::min(0.5,(t0 + bindingEnergy)/energy);
  G4double x2 = std::min(0.5,(tm + bindingEnergy)/energy);
  if(x1 >= x2) return tDelta;

  if(verbose > 1) {
    G4cout << "G4RDeIonisationSpectrum::SampleEnergy: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << G4endl;
  }

  // Access parameters
  G4DataVector p;

  // Access parameters
  for (G4int i=0; i<iMax; i++) 
  {
    G4double x = theParam->Parameter(Z, shell, i, e);
    if(i<4) x /= energy; 
    p.push_back(x);
  }

  if(p[3] > 0.5) p[3] = 0.5;

  G4double g = energy/electron_mass_c2 + 1.;
  p.push_back((2.0*g - 1.0)/(g*g));

  p[iMax-1] = Function(p[3], p);

  G4double aria1 = 0.0;
  G4double a1 = std::max(x1,p[1]);
  G4double a2 = std::min(x2,p[3]);
  if(a1 < a2) aria1 = IntSpectrum(a1, a2, p);
  G4double aria2 = 0.0;
  G4double a3 = std::max(x1,p[3]);
  G4double a4 = x2;
  if(a3 < a4) aria2 = IntSpectrum(a3, a4, p);

  G4double aria = (aria1 + aria2)*G4UniformRand();
  G4double amaj, fun, q, x, z1, z2, dx, dx1;

  //======= First aria to sample =====

  if(aria <= aria1) { 

    amaj = p[4];
    for (G4int j=5; j<iMax; j++) {
      if(p[j] > amaj) amaj = p[j];
    }

    a1 = 1./a1;
    a2 = 1./a2;

    G4int i;
    do {

      x = 1./(a2 + G4UniformRand()*(a1 - a2));
      z1 = p[1];
      z2 = p[3];
      dx = (p[2] - p[1]) / 3.0;
      dx1= std::exp(std::log(p[3]/p[2]) / 16.0);
      for (i=4; i<iMax-1; i++) {

        if (i < 7) {
          z2 = z1 + dx;
        } else if(iMax-2 == i) {
          z2 = p[3];
          break;
        } else {
          z2 = z1*dx1;
        }
        if(x >= z1 && x <= z2) break;
        z1 = z2;
      }
      fun = p[i] + (x - z1) * (p[i+1] - p[i])/(z2 - z1);

      if(fun > amaj) {
          G4cout << "WARNING in G4RDeIonisationSpectrum::SampleEnergy:" 
                 << " Majoranta " << amaj 
                 << " < " << fun
                 << " in the first aria at x= " << x
                 << G4endl;
      }

      q = amaj*G4UniformRand();

    } while (q >= fun);

  //======= Second aria to sample =====

  } else {

    amaj = std::max(p[iMax-1], Function(0.5, p)) * factor;
    a1 = 1./a3;
    a2 = 1./a4;

    do {

      x = 1./(a2 + G4UniformRand()*(a1 - a2));
      fun = Function(x, p);

      if(fun > amaj) {
          G4cout << "WARNING in G4RDeIonisationSpectrum::SampleEnergy:" 
                 << " Majoranta " << amaj 
                 << " < " << fun
                 << " in the second aria at x= " << x
                 << G4endl;
      }

      q = amaj*G4UniformRand();

    } while (q >= fun);

  }

  p.clear();

  tDelta = x*energy - bindingEnergy;

  if(verbose > 1) {
    G4cout << "tcut(MeV)= " << tMin/MeV 
           << "; tMax(MeV)= " << tMax/MeV 
           << "; x1= " << x1 
           << "; x2= " << x2 
           << "; a1= " << a1 
           << "; a2= " << a2 
           << "; x= " << x 
           << "; be= " << bindingEnergy 
           << "; e= " << e 
           << "; tDelta= " << tDelta 
           << G4endl;
  }


  return tDelta; 
}


G4double G4RDeIonisationSpectrum::IntSpectrum(G4double xMin, 
					    G4double xMax,
				      const G4DataVector& p) const
{
  // Please comment what IntSpectrum does
  G4double sum = 0.0;
  if(xMin >= xMax) return sum;

  G4double x1, x2, xs1, xs2, y1, y2, ys1, ys2, q;

  // Integral over interpolation aria
  if(xMin < p[3]) {

    x1 = p[1];
    y1 = p[4];

    G4double dx = (p[2] - p[1]) / 3.0;
    G4double dx1= std::exp(std::log(p[3]/p[2]) / 16.0);

    for (size_t i=0; i<19; i++) {

      q = 0.0;
      if (i < 3) {
        x2 = x1 + dx;
      } else if(18 == i) {
        x2 = p[3];
      } else {
        x2 = x1*dx1;
      }

      y2 = p[5 + i];

      if (xMax <= x1) {
        break;
      } else if (xMin < x2) {

        xs1 = x1;
        xs2 = x2;
        ys1 = y1;
        ys2 = y2;

        if (x2 > x1) {
          if (xMin > x1) {
            xs1 = xMin;
            ys1 += (xs1 - x1)*(y2 - y1)/(x2 - x1);
	  } 
          if (xMax < x2) {
            xs2 = xMax;
            ys2 += (xs2 - x2)*(y1 - y2)/(x1 - x2);
	  } 
          if (xs2 > xs1) {
            q = (ys1*xs2 - ys2*xs1)/(xs1*xs2) 
              +  std::log(xs2/xs1)*(ys2 - ys1)/(xs2 - xs1);
            sum += q;
            if(p.size() == 26) G4cout << "i= " << i << "  q= " << q << " sum= " << sum << G4endl;
	  }
	}  
      }
      x1 = x2;
      y1 = y2;
    }
  }

  // Integral over aria with parametrised formula 

  x1 = std::max(xMin, p[3]);
  if(x1 >= xMax) return sum;
  x2 = xMax;

  xs1 = 1./x1;
  xs2 = 1./x2;
  q = (xs1 - xs2)*(1.0 - p[0]) 
       - p[iMax]*std::log(x2/x1) 
       + (1. - p[iMax])*(x2 - x1)
       + 1./(1. - x2) - 1./(1. - x1) 
       + p[iMax]*std::log((1. - x2)/(1. - x1))
       + 0.25*p[0]*(xs1*xs1 - xs2*xs2);
  sum += q;
  if(p.size() == 26) G4cout << "param...  q= " << q <<  " sum= " << sum << G4endl;

  return sum;
} 


G4double G4RDeIonisationSpectrum::AverageValue(G4double xMin, 
				             G4double xMax,
				       const G4DataVector& p) const
{
  G4double sum = 0.0;
  if(xMin >= xMax) return sum;

  G4double x1, x2, xs1, xs2, y1, y2, ys1, ys2;

  // Integral over interpolation aria
  if(xMin < p[3]) {

    x1 = p[1];
    y1 = p[4];

    G4double dx = (p[2] - p[1]) / 3.0;
    G4double dx1= std::exp(std::log(p[3]/p[2]) / 16.0);

    for (size_t i=0; i<19; i++) {

      if (i < 3) {
        x2 = x1 + dx;
      } else if(18 == i) {
        x2 = p[3];
      } else {
        x2 = x1*dx1;
      }

      y2 = p[5 + i];

      if (xMax <= x1) {
        break;
      } else if (xMin < x2) {

        xs1 = x1;
        xs2 = x2;
        ys1 = y1;
        ys2 = y2;

        if (x2 > x1) {
          if (xMin > x1) {
            xs1 = xMin;
            ys1 += (xs1 - x1)*(y2 - y1)/(x2 - x1);
	  } 
          if (xMax < x2) {
            xs2 = xMax;
            ys2 += (xs2 - x2)*(y1 - y2)/(x1 - x2);
	  } 
          if (xs2 > xs1) {
            sum += std::log(xs2/xs1)*(ys1*xs2 - ys2*xs1)/(xs2 - xs1) 
                +  ys2 - ys1;
	  }
	}  
      }
      x1 = x2;
      y1 = y2;

    }
  }

  // Integral over aria with parametrised formula 

  x1 = std::max(xMin, p[3]);
  if(x1 >= xMax) return sum;
  x2 = xMax;

  xs1 = 1./x1;
  xs2 = 1./x2;

  sum  += std::log(x2/x1)*(1.0 - p[0]) 
        + 0.5*(1. - p[iMax])*(x2*x2 - x1*x1)
        + 1./(1. - x2) - 1./(1. - x1) 
        + (1. + p[iMax])*std::log((1. - x2)/(1. - x1))
        + 0.5*p[0]*(xs1 - xs2);

  return sum;
} 


void G4RDeIonisationSpectrum::PrintData() const 
{
  theParam->PrintData();
}

G4double G4RDeIonisationSpectrum::MaxEnergyOfSecondaries(G4double kineticEnergy,
						       G4int, // Z = 0,
						       const G4ParticleDefinition* ) const
{
  return 0.5 * kineticEnergy;
}
