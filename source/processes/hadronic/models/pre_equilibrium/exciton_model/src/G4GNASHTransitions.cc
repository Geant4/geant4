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
// $Id$
//
// 20.08.2010 V.Ivanchenko move constructor and destructor to the source 

#include "G4GNASHTransitions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PreCompoundParameters.hh"
#include "G4HadronicException.hh"

G4GNASHTransitions::G4GNASHTransitions()
{}

G4GNASHTransitions::~G4GNASHTransitions()
{}

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
  if ( x < 2.0*MeV ) theMatrixElement *= x/std::sqrt(14.0*MeV*MeV);
  else if ( x < 7.0*MeV ) x *= std::sqrt(x/7.0*MeV);
  else if ( x < 15.0*MeV ) ;
  else x *= std::sqrt(15.0*MeV/x);

  // gg = (6.0/pi2)*a*A
  G4double gg =  (6.0/pi2)*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*A;

  G4double Epauli = ((P+1.0)*(P+1.0) + (H+1.0)*(H+1.0) + (P+1.0) - 3.0*(H-1.0))/4.0;

  G4double Probability = gg*gg*gg *(E-Epauli)*(E-Epauli);
  Probability /= 2.0*(N+1.0)*h_Planck;
  Probability *= theMatrixElement;

  return Probability;
}

void G4GNASHTransitions::PerformTransition(G4Fragment & result)
{
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
}
