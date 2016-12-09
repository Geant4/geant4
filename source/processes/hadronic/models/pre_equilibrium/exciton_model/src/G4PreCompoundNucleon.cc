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
// $Id: G4PreCompoundNucleon.cc 100378 2016-10-19 15:03:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundNucleon
//
// Author:         V.Lara
//
// Modified:  
// 10.02.2009 J. M. Quesada cleanup
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//

#include "G4PreCompoundNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4PreCompoundNucleon::G4PreCompoundNucleon(
    const G4ParticleDefinition* part, G4VCoulombBarrier* aCoulombBarrier)
  : G4PreCompoundFragment(part,aCoulombBarrier)
{}

G4PreCompoundNucleon::~G4PreCompoundNucleon()
{}

G4double G4PreCompoundNucleon::
ProbabilityDistributionFunction(G4double eKin, 
				const G4Fragment& aFragment)
{
  G4double U = aFragment.GetExcitationEnergy();
  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();
  G4int N = P + H;

  static const G4double sixoverpi2 = 6.0/CLHEP::pi2;
  G4double g0 = sixoverpi2*theFragA*theParameters->GetLevelDensity();
  G4double g1 = sixoverpi2*theResA*theParameters->GetLevelDensity();
  
  G4double A0 = G4double(P*P+H*H+P-3*H)/(4.0*g0);
  G4double A1 = (A0 - 0.5*P)/g1;  

  G4double E0 = U - A0;
  if (E0 <= 0.0) { return 0.0; }

  G4double E1 = U - eKin -  theBindingEnergy - A1;
  if (E1 <= 0.0) { return 0.0; }

  G4double rj = GetRj(P, aFragment.GetNumberOfCharged());
  G4double xs = CrossSection(eKin);

  if (rj <0.0 || xs < 0.0) { return 0.0; }

  static const G4double fact = 2*CLHEP::millibarn
    /(CLHEP::pi2*CLHEP::hbarc*CLHEP::hbarc*CLHEP::hbarc);
  G4double Probability = fact * theReducedMass * rj * xs * eKin * P * (N-1) 
    * g4calc->powN(g1*E1/(g0*E0),N-2) * g1 / (E0*g0*g0);
  
  return Probability;
}
















