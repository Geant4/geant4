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
// $Id: G4eIonisationSpectrum.cc,v 1.7 2001-11-08 10:39:11 vnivanch Exp $
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
//
// -------------------------------------------------------------------
//

#include "G4eIonisationSpectrum.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4eIonisationParameters.hh"
#include "G4DataVector.hh"
#include "Randomize.hh"


G4eIonisationSpectrum::G4eIonisationSpectrum(): 
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

  // Access parameters
  G4double t1 = theParam->Parameter(Z, shell, 14, 0.0);
  G4double t2 = theParam->Parameter(Z, shell, 15, 0.0);

  if(verbose > 1) {
    G4cout << "G4eIonisationSpectrum::Probability: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; t1= " << t1 
           << "; t2= " << t2 
           << G4endl;
  }

  G4int i;
  G4int iMax = 14;
  G4DataVector p;

  for (i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }

  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();
    
  G4double val = 0.0;
  G4double nor = 0.0;
  
  // Integration is performed in 3 energy arias, in each probability 
  // is integrated between a2 and a3, normalisation - between a1 and a3

  // --- First function ---

  G4double a1 = lowestE; 
  G4double a2 = G4std::min(t0,t1); 
  G4double a3 = G4std::min(eMax,t1);
  G4double a4 = G4std::min(tm,t1);

  if(a2 < a4) val += IntSpectrum(6, a2, a4, bindingEnergy, p);
  if(a1 < a3) nor += IntSpectrum(6, a1, a3, bindingEnergy, p);

  // --- Second function ---

  if(eMax > t1 && t2 > t1) {
    a1 = t1; 
    a3 = G4std::min(eMax, t2);
    a2 = G4std::max(t0,t1); 
    a4 = G4std::min(a3,G4std::max(a2,tm)); 
    G4double c1 = theParam->Parameter(Z, shell, 12, e);
    G4double c2 = theParam->Parameter(Z, shell, 13, e);

    if(a2 < a4) val += c1 * (pow(a4, c2 + 1.) - pow(a2, c2 + 1.)) / (c2 + 1.);
    if(a1 < a3) nor += c1 * (pow(a3, c2 + 1.) - pow(a1, c2 + 1.)) / (c2 + 1.);
  }

  // --- Third function ---

  if(eMax > t2) {

    a1 = t2; 
    a3 = eMax;
    a2 = G4std::max(t0,t2); 
    a4 = G4std::min(a3,G4std::max(a2,tm)); 

    if(a2 < a3) val += IntSpectrum(4, a2, a3, bindingEnergy, p);
    if(a1 < a3) nor += IntSpectrum(4, a1, a3, bindingEnergy, p);
  }
  if(verbose > 1) {
    G4cout << "tcut= " << tMin 
           << "; tMax= " << tMax 
           << "; t0= " << t0 
           << "; tm= " << tm 
           << "; eMax= " << eMax 
           << "; t1= " << t1 
           << "; t2= " << t2 
           << "; val= " << val 
           << "; nor= " << nor 
           << G4endl;
  }
  p.clear();

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

  // Access parameters
  G4double t1 = theParam->Parameter(Z, shell, 14, 0.0);
  G4double t2 = theParam->Parameter(Z, shell, 15, 0.0);

  if(verbose > 1) {
    G4cout << "G4eIonisationSpectrum::AverageEnergy: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << "; t1= " << t1 
           << "; t2= " << t2 
           << G4endl;
  }

  G4int i;
  G4int iMax = 14;
  G4DataVector p;

  for (i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }
 
  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
    Shell(Z, shell)->BindingEnergy();
    
  G4double val = 0.0;
  G4double nor = 0.0;
  
  // Integration is performed in 3 energy arias, in each probability 
  // is integrated between a2 and a3, normalisation - between a1 and a3

  // --- First function ---

  G4double a1 = lowestE; 
  G4double a2 = G4std::min(t0,t1); 
  G4double a3 = G4std::min(eMax,t1);
  G4double a4 = G4std::min(tm,t1);

  if(a2 < a4) val += AverageValue(6, a2, a4, bindingEnergy, p);
  if(a1 < a3) nor += IntSpectrum (6, a1, a3, bindingEnergy, p);

  // --- Second function ---

  if(eMax > t1 && t2 > t1) {
    a1 = t1; 
    a3 = G4std::min(eMax, t2);
    a2 = G4std::min(G4std::max(t0,t1),eMax); 
    a4 = G4std::min(G4std::max(t0,t1),tm); 
    G4double c1 = theParam->Parameter(Z, shell, 12, e);
    G4double c2 = theParam->Parameter(Z, shell, 13, e);

    if(a2 < a4) {
      val += c1 * (pow(a4, c2 + 2.) - pow(a2, c2 + 2.)) / (c2 + 2.);
      val += c1 * bindingEnergy *
                  (pow(a4, c2 + 1.) - pow(a2, c2 + 1.)) / (c2 + 1.);
    }
    if(a1 < a3) nor += c1 * (pow(a3, c2 + 1.) - pow(a1, c2 + 1.)) / (c2 + 1.);
  }

  // --- Third function ---

  if(eMax > t2) {

    a1 = t2; 
    a3 = eMax;
    a2 = G4std::max(t0,t2); 
    a4 = G4std::min(a3,G4std::max(a2,tm)); 

    if(a2 < a4) val += AverageValue(4, a2, a4, bindingEnergy, p);
    if(a1 < a3) nor += IntSpectrum (4, a1, a3, bindingEnergy, p);
  }
  if(verbose > 1) {
    G4cout << "tcut= " << tMin 
           << "; tMax= " << tMax 
           << "; t0= " << t0 
           << "; tm= " << tm 
           << "; eMax= " << eMax 
           << "; t1= " << t1 
           << "; t2= " << t2 
           << "; val= " << val 
           << "; nor= " << nor 
           << G4endl;
  }
  p.clear();

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

  if(verbose > 1 && Z == 14) {
    G4cout << "G4eIonisationSpectrum::Probability: Z= " << Z
           << "; shell= " << shell
           << "; E(keV)= " << e/keV
           << G4endl;
  }

  // Access parameters
  G4double t1 = theParam->Parameter(Z, shell, 14, 0.0);
  G4double t2 = theParam->Parameter(Z, shell, 15, 0.0);
  G4int i;
  G4int iMax = 14;
  G4DataVector p;
  for (i=0; i<iMax; i++) 
  {
    p.push_back(theParam->Parameter(Z, shell, i, e));
  }
  G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4double bindingEnergy = (transitionManager->Shell(Z, shell))->BindingEnergy();

  G4double aria1 = 0.0;
  G4double a1 = G4std::min(t0,t1);
  G4double a2 = G4std::min(tm,t1);
  if(a1 < a2) aria1 = IntSpectrum(6, a1, a2, bindingEnergy, p);
  G4double aria2 = 0.0;
  a1 = G4std::max(t0,t1);
  a2 = G4std::min(tm,t2);
  G4double c2 = p[13] + 1.;
  if(a1 < a2) aria2 = p[12] * (pow(a2, c2) - pow(a1, c2)) /c2;
  G4double aria3 = 0.0;
  a1 = G4std::max(t0,t2);
  a2 = G4std::max(tm,t2);
  if(a1 < a2) aria3 = IntSpectrum(4, a1, a2, bindingEnergy, p);

  G4double aria = (aria1 + aria2 + aria3)*G4UniformRand();
  G4double amaj, aleft, fun, q, q1, q2, q3;

  //======= First function to sample =====

  if(aria <= aria1) { 

    a1 = G4std::min(t0,t1);
    a2 = G4std::min(tm,t1);
    amaj = MaxFunction(6, a1, a2, bindingEnergy, p);
    aleft = 1.21*Function(6, a1, bindingEnergy, p);

    if(aleft > amaj) {

      aleft *= 2.0;
      q3 = (a1 + bindingEnergy);
      aleft *= q3*q3; 
      q1 = 1./q3;
      q2 = 1./(a2 + bindingEnergy);

      do {

        q3 = q1 + G4UniformRand()*(q2 - q1);
        tDelta = 1./q3 - bindingEnergy;
        fun  = Function(6, tDelta, bindingEnergy, p)/(q3*q3);

        if(fun > aleft) {
          G4cout << "WARNING in G4eIonisationSpectrum::SampleEnergy:" 
                 << " 1st majoranta " << aleft 
                 << " < " << fun
                 << G4endl;
        }

        q = aleft*G4UniformRand();

      } while (q >= fun);

    } else {

      do {

        tDelta = a1 + G4UniformRand()*(a2 - a1);
        fun  = Function(6, tDelta, bindingEnergy, p);

        if(fun > amaj) {
          G4cout << "WARNING in G4eIonisationSpectrum::SampleEnergy:" 
                 << " 1-st majoranta " << amaj 
                 << " < " << fun
                 << G4endl;
        }

        q = amaj*G4UniformRand();

      } while (q >= fun);
    }

  //======= Second function to sample =====

  } else if(aria < aria1 + aria2) { 

    a1 = G4std::max(t0,t1);
    a2 = G4std::min(tm,t2);
    a1 = pow(a1, c2);
    a2 = pow(a2, c2);
    q  = a1 + G4UniformRand()*(a2 - a1);
    tDelta = pow(q, 1.0/c2);

  //======= Third function to sample =====

   } else { 

    a1 = G4std::max(t0,t2);
    a2 = G4std::max(tm,t2);
    amaj = MaxFunction(4, a1, a2, bindingEnergy, p);
    aleft = 1.21*Function(4, a1, bindingEnergy, p);

    if(aleft > amaj) {

      aleft *= 2.0;
      q3 = (a1 + bindingEnergy);
      aleft *= q3*q3; 
      q1 = 1./q3;
      q2 = 1./(a2 + bindingEnergy);

      do {

        q3 = q1 + G4UniformRand()*(q2 - q1);
        tDelta = 1./q3 - bindingEnergy;
        fun  = Function(4, tDelta, bindingEnergy, p)/(q3*q3);

        if(fun > aleft) {
          G4cout << "WARNING in G4eIonisationSpectrum::SampleEnergy:" 
                 << " 3st majoranta " << aleft 
                 << " < " << fun
                 << G4endl;
        }

        q = aleft*G4UniformRand();

      } while (q >= fun);

    } else {

      do {

        tDelta = a1 + G4UniformRand()*(a2 - a1);
        fun  = Function(4, tDelta, bindingEnergy, p);

        if(fun > amaj) {
          G4cout << "WARNING in G4eIonisationSpectrum::SampleEnergy:" 
                 << " 3-d majoranta " << amaj 
                 << " < " << fun
                 << G4endl;
        }

        q = amaj*G4UniformRand();

      } while (q >= fun);
    }
  }

  p.clear();

  return tDelta; 
}


G4double G4eIonisationSpectrum::IntSpectrum(size_t n, 
					    G4double tMin, 
					    G4double tMax,
					    G4double b,
					    const G4DataVector& p) const
{
  // Please comment what IntSpectrum does
  G4int k = 0;
  if(n == 4) k = 7;
  G4double zMax = 1./(tMax + b);
  G4double zMin = 1./(tMin + b);
  G4double yMax = zMax;
  G4double yMin = zMin;
  G4double x = p[k]*(yMin - yMax);
  for(size_t i=1; i<n; i++) {
    yMax *= zMax;
    yMin *= zMin;
    x += p[i + k]*(yMin - yMax)/(1.0 + (G4double)i);
  }
  if(x < 0.0) x = 0.0;
  return x;
} 


G4double G4eIonisationSpectrum::AverageValue(size_t n, 
					     G4double tMin, 
					     G4double tMax,
					     G4double b,
					     const G4DataVector& p) const
{
  G4int k = 0;
  if(n == 4) k = 7;
  G4double zMax = 1./(tMax + b);
  G4double zMin = 1./(tMin + b);
  G4double yMax = zMax;
  G4double yMin = zMin;
  G4double x = p[k]*log(zMin/zMax) + p[k + 1]*(yMin - yMax);
  for(size_t i=2; i<n; i++) {
    yMax *= zMax;
    yMin *= zMin;
    x += p[k + i]*(yMin - yMax)/((G4double)i);
  }
  if(x < 0.0) x = 0.0;
  return x;
} 


G4double G4eIonisationSpectrum::Function(size_t n, 
					 G4double e, 
					 G4double b,
					 const G4DataVector& p) const
{
  // Please comment what Function does
  G4int k = 0;
  if(n == 4) k = 7;
  G4double z = 1./(e + b);
  G4double y = z;
  G4double x = 0.0;
  for(size_t i=0; i<n; i++) {
    y *= z;
    x += p[k + i]*y;
  }
  if(x < 0.0) x = 0.0;
  return x;
} 


G4double G4eIonisationSpectrum::MaxFunction(size_t n, 
					    G4double tMin, 
					    G4double tMax,
					    G4double b,
					    const G4DataVector& p) const
{
  // Please comment what MaxFunction does

  size_t nbin  = 11;
  G4double de = 0.1*(tMax - tMin);
  G4double x = 0.0;
  G4double y, e;
  for(size_t i=0; i<nbin; i++) {
    e = tMin + de * ((G4double)i);
    y = Function(n, e, b, p);
    if(y > x) x = y;
  }
  if(x < 0.0) x = 0.0;
  return x*1.2;
} 

void G4eIonisationSpectrum::PrintData() const 
{
  theParam->PrintData();
}





