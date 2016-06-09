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
//
// $Id: G4VPreCompoundIon.cc,v 1.15 2002/12/12 19:17:33 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// by V. Lara

#include "G4VPreCompoundIon.hh"
#include "G4PreCompoundParameters.hh"
//#include "G4EvaporationLevelDensityParameter.hh"


G4double G4VPreCompoundIon::
ProbabilityDistributionFunction(const G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) return 0.0;

  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();

  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;

  G4double A0 = (P*P+H*H+P-H)/4.0 - H/2.0;
  G4double A1 = A0 + GetA()*(GetA()-2.0*P-1.0)/4.0;
  G4double Aj = GetA()*(GetA()+1.0)/4.0;

  G4double E0 = G4std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  G4double E1 = G4std::max(0.0,U - eKin - GetBindingEnergy() - A1);
  G4double Ej = G4std::max(0.0,U + GetBindingEnergy() - Aj);

  // g = 0.595*a*A
//  G4EvaporationLevelDensityParameter theLDP;
  G4double g0 = 0.595*aFragment.GetA() * 
      G4PreCompoundParameters::GetAddress()->GetLevelDensity();
//      theLDP.LevelDensityParameter(aFragment.GetA(),aFragment.GetZ(),U);
  G4double g1 = 0.595*GetRestA() * 
      G4PreCompoundParameters::GetAddress()->GetLevelDensity();
//      theLDP.LevelDensityParameter(GetRestA(),GetRestZ(),U);
  G4double gj = 0.595*GetA() *
      G4PreCompoundParameters::GetAddress()->GetLevelDensity();
//      theLDP.LevelDensityParameter(GetA(),GetZ(),U);

  G4double pA = (3.0/4.0) * sqrt(G4std::max(0.0, 2.0/(GetReducedMass()*(eKin+GetBindingEnergy()))))*
      GetAlpha()*(eKin+GetBeta())/(r0*pow(GetRestA(),1.0/3.0)) *
      CoalescenceFactor(aFragment.GetA()) * FactorialFactor(N,P);

  G4double pB = pow((g1*E1)/(g0*E0),N-GetA()-1.0)*(g1/g0);

  G4double pC = pow((gj*Ej)/(g0*E0),GetA()-1.0)*(gj/g0)/E0;

  G4double Probability = pA * pB * pC;

  return Probability;
}







