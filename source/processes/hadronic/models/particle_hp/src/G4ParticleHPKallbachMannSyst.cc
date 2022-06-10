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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080801 Protect div0 error, when theCompundFraction is 1 by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
// June-2019 - E. Mendoza --> perform some corrections

#include "G4ParticleHPKallbachMannSyst.hh" 
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" 
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4HadronicException.hh" 

G4double G4ParticleHPKallbachMannSyst::Sample(G4double anEnergy)
{
  G4double result=0.;
  
  G4double zero = GetKallbachZero(anEnergy);
  if(zero>1) zero=1.;
  if(zero<-1)zero=-1.;
  G4double max = Kallbach(zero, anEnergy);
  G4double upper = Kallbach(1., anEnergy);
  G4double lower = Kallbach(-1., anEnergy);
  if(upper>max) max=upper;
  if(lower>max) max=lower;
  G4double value, random;

  G4int icounter=0;
  G4int icounter_max=1024;
  do
  {
    icounter++;
    if ( icounter > icounter_max ) {
       G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
       break;
    }
    result = 2.*G4UniformRand()-1;  
    value = Kallbach(result, anEnergy)/max;
    random = G4UniformRand();
  }
  while(random>value); // Loop checking, 11.05.2015, T. Koi
  
  return result;
}

G4double G4ParticleHPKallbachMannSyst::Kallbach(G4double cosTh, G4double anEnergy)
{
  // Kallbach-Mann systematics without normalization.
  G4double result;
  G4double theX = A(anEnergy)*cosTh;
  result = 0.5*(G4Exp( theX)*(1+theCompoundFraction)
               +G4Exp(-theX)*(1-theCompoundFraction));
  return result;
}

G4double G4ParticleHPKallbachMannSyst::GetKallbachZero(G4double anEnergy)
{
  G4double result;
  //delta 2.0e-16 in not good.
  //delta 4.0e-16 is OK 
  //safety factor of 2
  G4double delta = 8.0e-16;
  if ( std::abs (theCompoundFraction - 1 ) < delta ) { 
     theCompoundFraction = 1.0-delta;   
  } 
  result = 0.5 * (1./A(anEnergy)) * G4Log((1-theCompoundFraction)/(1+theCompoundFraction));
  return result;
}

G4double G4ParticleHPKallbachMannSyst::A(G4double anEnergy)
{ 
  G4double result;
  G4double C1 = 0.04/MeV;
  G4double C2 = 1.8E-6/(MeV*MeV*MeV);
  G4double C3 = 6.7E-7/(MeV*MeV*MeV*MeV);
  
  G4double epsa = anEnergy*theTargetMass/(theTargetMass+theIncidentMass);
  G4int Ac = theTargetA+theProjectileA; 
  G4int Nc = Ac - theTargetZ-theProjectileZ;
  G4int AA = theTargetA;
  G4int ZA = theTargetZ;
  G4double ea = epsa+SeparationEnergy(Ac, Nc, AA, ZA,theProjectileA,theProjectileZ);
  G4double Et1 = 130*MeV;
  G4double R1 = std::min(ea, Et1);
  // theProductEnergy is still in CMS!!!
  G4double epsb = theProductEnergy*(theProductMass+theResidualMass)/theResidualMass;
  G4int AB = theResidualA;
  G4int ZB = theResidualZ;
  G4double eb = epsb+SeparationEnergy(Ac, Nc, AB, ZB,theProductA, theProductZ);
  G4double X1 = R1*eb/ea;
  G4double Et3 = 41*MeV;
  G4double R3 = std::min(ea, Et3);
  G4double X3 = R3*eb/ea;

  G4double Ma=1;
  G4double mb=1;
  if(theProjectileA==1 || (theProjectileZ==1 && theProjectileA==2)){Ma=1;}//neutron,proton,deuteron
  else if(theProjectileA==4 && theProjectileZ==2){Ma=0;}//alpha
  else if(theProjectileA==3 && (theProjectileZ==1 || theProjectileZ==2)){Ma=0.5;}//tritum,He3 : set intermediate value
  else
  {
    throw G4HadronicException(__FILE__, __LINE__, "Severe error in the sampling of Kallbach-Mann Systematics");
  }
  if(theProductA==1 && theProductZ==0){mb=1./2.;}//neutron
  else if(theProductA==4 && theProductZ==2){mb=2;}//alpha
  else{mb=1;}

  result = C1*X1 + C2*G4Pow::GetInstance()->powN(X1, 3) + C3*Ma*mb*G4Pow::GetInstance()->powN(X3, 4);
  return result;

}

G4double G4ParticleHPKallbachMannSyst::SeparationEnergy(G4int Ac, G4int Nc, G4int AA, G4int ZA,G4int Abinding,G4int Zbinding)
{
  G4double result;
  G4int NA = AA-ZA;
  G4int Zc = Ac-Nc;
  result = 15.68*(Ac-AA);
  result += -28.07*((Nc-Zc)*(Nc-Zc)/(G4double)Ac - (NA-ZA)*(NA-ZA)/(G4double)AA);
  result += -18.56*(G4Pow::GetInstance()->A23(G4double(Ac)) - G4Pow::GetInstance()->A23(G4double(AA)));
  result +=  33.22*((Nc-Zc)*(Nc-Zc)/G4Pow::GetInstance()->powA(G4double(Ac), 4./3.) - (NA-ZA)*(NA-ZA)/G4Pow::GetInstance()->powA(G4double(AA), 4./3.));
  result += -0.717*(Zc*Zc/G4Pow::GetInstance()->A13(G4double(Ac))-ZA*ZA/G4Pow::GetInstance()->A13(G4double(AA)));
  result +=  1.211*(Zc*Zc/(G4double)Ac-ZA*ZA/(G4double)AA);
  G4double totalBinding(0);
  if(Zbinding==0&&Abinding==1) totalBinding=0;
  if(Zbinding==1&&Abinding==1) totalBinding=0;
  if(Zbinding==1&&Abinding==2) totalBinding=2.224596;
  if(Zbinding==1&&Abinding==3) totalBinding=8.481798;
  if(Zbinding==2&&Abinding==3) totalBinding=7.718043;
  if(Zbinding==2&&Abinding==4) totalBinding=28.29566;
  result += -totalBinding;
  result *= MeV;
  return result;
}
