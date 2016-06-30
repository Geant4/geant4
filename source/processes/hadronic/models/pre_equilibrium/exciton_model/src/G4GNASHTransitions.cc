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
// $Id: G4GNASHTransitions.cc 96603 2016-04-25 13:29:51Z gcosmo $
//
// 20.08.2010 V.Ivanchenko move constructor and destructor to the source 

#include "G4GNASHTransitions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4HadronicException.hh"
#include "Randomize.hh"

G4GNASHTransitions::G4GNASHTransitions()
{
  theParameters = G4NuclearLevelData::GetInstance()->GetParameters();
}

G4GNASHTransitions::~G4GNASHTransitions()
{}

G4double G4GNASHTransitions::
CalculateProbability(const G4Fragment & aFragment)
{
  static const G4double k = 135.0 *CLHEP::MeV*CLHEP::MeV*CLHEP::MeV;
  G4double E = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  G4double A = aFragment.GetA_asInt();

  G4double theMatrixElement(k*N/(A*A*A*E));
  G4double x = E/(N*CLHEP::MeV);
  static const G4double xf = std::sqrt(2.0/7.0);
  if ( x < 2.0)      { x *= xf;   }
  else if ( x < 7.0) { x *= std::sqrt(x/7.0);  }
  else if ( x > 15.0){ x *= std::sqrt(15.0/x); }
  theMatrixElement *= x;

  G4double gg =  (6.0/pi2)*theParameters->GetLevelDensity()*A;

  G4double Epauli = ((P+1.0)*(P+1.0) + (H+1.0)*(H+1.0) + (P+1.0) - 3.0*(H-1.0))*0.25;

  G4double Probability = gg*gg*gg *(E-Epauli)*(E-Epauli);
  Probability *= theMatrixElement/(2.0*(N+1.0)*CLHEP::h_Planck);

  return Probability;
}

void G4GNASHTransitions::PerformTransition(G4Fragment & result)
{
  result.SetNumberOfParticles(result.GetNumberOfParticles()+1);
  result.SetNumberOfHoles(result.GetNumberOfHoles()+1);
  if (G4UniformRand()*result.GetA_asInt() <= G4double(result.GetZ_asInt()))
    {
      result.SetNumberOfCharged(result.GetNumberOfCharged()+1);
    }

  if (result.GetNumberOfParticles() < result.GetNumberOfCharged())
    {
      result.SetNumberOfCharged(result.GetNumberOfParticles());
    }
}
