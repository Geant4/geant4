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
// $Id: G4VPreCompoundNucleon.cc,v 1.5 2003/06/16 17:07:34 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// by V. Lara

#include "G4VPreCompoundNucleon.hh"
#include "G4PreCompoundParameters.hh"
//#include "G4EvaporationLevelDensityParameter.hh"


G4double G4VPreCompoundNucleon::
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
  G4double A1 = A0 - P/2.0;

  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  G4double E1 = std::max(0.0,U - eKin - GetBindingEnergy() - A1);

  // g = (6.0/pi2)*a*A
  //  G4EvaporationLevelDensityParameter theLDP;
  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(G4int(aFragment.GetA()),G4int(aFragment.GetZ()),U);
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
    //    theLDP.LevelDensityParameter(G4int(GetRestA()),G4int(GetRestZ()),U);


  G4double Probability = 2.0/(hbarc*hbarc*hbarc) * GetReducedMass() * 
      r0 * r0 * pow(GetRestA(),2.0/3.0) * GetAlpha() * (eKin + GetBeta()) *
      P*(N-1.0) * pow(g1*E1/(g0*E0),N-2.0)/E0 *
      g1/(g0*g0);

  return Probability;
}


