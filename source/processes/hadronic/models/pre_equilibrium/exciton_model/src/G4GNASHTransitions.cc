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
#include "G4GNASHTransitions.hh"
#include "G4PreCompoundParameters.hh"
#include "G4HadronicException.hh"

const G4GNASHTransitions & G4GNASHTransitions::
operator=(const G4GNASHTransitions & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GNASHTransitions::operator= meant to not be accesable");
  return *this;
}

G4bool G4GNASHTransitions::operator==(const G4GNASHTransitions & ) const
{
  return false;
}

G4bool G4GNASHTransitions::operator!=(const G4GNASHTransitions & ) const
{
  return true;
}


G4double G4GNASHTransitions::
CalculateProbability(const G4Fragment & aFragment)
{
  const G4double k = 135.0 * MeV*MeV*MeV;
  G4double E = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  G4double A = aFragment.GetA();

  G4double theMatrixElement(k*N/(A*A*A*E));
  G4double x = E/N;
  if ( x < 2.0*MeV ) theMatrixElement *= x/sqrt(14.0*MeV*MeV);
  else if ( x < 7.0*MeV ) x *= sqrt(x/7.0*MeV);
  else if ( x < 15.0*MeV ) ;
  else x *= sqrt(15.0*MeV/x);

  // g = (6.0/pi2)*a*A
  G4double g =  (6.0/pi2)*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*A;

  G4double Epauli = ((P+1.0)*(P+1.0) + (H+1.0)*(H+1.0) + (P+1.0) - 3.0*(H-1.0))/4.0;

  G4double Probability = g*g*g *(E-Epauli)*(E-Epauli);
  Probability /= 2.0*(N+1.0)*h_Planck;
  Probability *= theMatrixElement;


  return Probability;
}

G4Fragment G4GNASHTransitions::
PerformTransition(const G4Fragment & aFragment)
{
  G4Fragment result(aFragment);
  result.SetNumberOfParticles(result.GetNumberOfParticles()+1);
  result.SetNumberOfHoles(result.GetNumberOfHoles()+1);
  if (G4UniformRand() <= result.GetZ()/result.GetA())
    {
      result.SetNumberOfCharged(result.GetNumberOfCharged()+1);
    }

  if (result.GetNumberOfParticles() < result.GetNumberOfCharged())
    {
      result.SetNumberOfCharged(result.GetNumberOfParticles());
    }

  return result;
}
