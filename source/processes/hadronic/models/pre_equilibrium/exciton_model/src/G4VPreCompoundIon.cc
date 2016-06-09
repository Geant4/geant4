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
// $Id: G4VPreCompoundIon.cc,v 1.2 2004/12/07 13:50:44 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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
  G4double A1 = std::max(0.0,A0 + GetA()*(GetA()-2.0*P-1.0)/4.0);
  G4double Aj = GetA()*(GetA()+1.0)/4.0;

  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  G4double E1 = std::max(0.0,GetMaximalKineticEnergy() + GetCoulombBarrier() - eKin - A1);
  G4double Ej = std::max(0.0,U + GetBindingEnergy() - Aj);

  // g = (6.0/pi2)*a*A
  //  G4EvaporationLevelDensityParameter theLDP;
  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(static_cast<G4int>(aFragment.GetA()),static_cast<G4int>(aFragment.GetZ()),U);
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(GetRestA(),GetRestZ(),U);
  G4double gj = (6.0/pi2)*GetA() *
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(GetA(),GetZ(),U);

  G4double pA = (3.0/4.0) * std::sqrt(std::max(0.0, 2.0/(GetReducedMass()*(eKin+GetBindingEnergy()))))*
      GetAlpha()*(eKin+GetBeta())/(r0*std::pow(GetRestA(),1.0/3.0)) *
      CoalescenceFactor(aFragment.GetA()) * FactorialFactor(N,P);

  G4double pB = std::pow((g1*E1)/(g0*E0),N-GetA()-1.0)*(g1/g0);

  G4double pC = std::pow((gj*Ej)/(g0*E0),GetA()-1.0)*(gj/g0)/E0;

  G4double Probability = pA * pB * pC;

  return Probability;
}







