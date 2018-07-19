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
// $Id: G4PreCompoundIon.cc 100378 2016-10-19 15:03:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundIon
//
// Author:         V.Lara
//
// Modified:  
// 10.02.2009 J. M. Quesada fixed bug in level density of light fragments  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//

#include "G4PreCompoundIon.hh"
#include "G4PhysicalConstants.hh"

G4PreCompoundIon::
G4PreCompoundIon(const G4ParticleDefinition* part,
		 G4VCoulombBarrier* aCoulombBarrier)
  : G4PreCompoundFragment(part,aCoulombBarrier)
{
  G4double r0 = theParameters->GetR0();
  fact = 0.75*CLHEP::millibarn/(CLHEP::pi*r0*r0*r0);
}

G4PreCompoundIon::~G4PreCompoundIon()
{}

G4double G4PreCompoundIon::
ProbabilityDistributionFunction(G4double eKin, 
				const G4Fragment& aFragment)
{
  G4double efinal = eKin + theBindingEnergy;
  if(efinal <= 0.0 ) { return 0.0; } 

  G4double U = aFragment.GetExcitationEnergy();
  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();
  G4int A = GetA();
  G4int N = P + H;

  static const G4double sixoverpi2 = 6.0/CLHEP::pi2;
  G4double g0 = sixoverpi2*theFragA*theParameters->GetLevelDensity();
  G4double g1 = sixoverpi2*theResA*theParameters->GetLevelDensity();

  G4double gj = g1;

  G4double A0 = G4double(P*P+H*H+P-3*H)/(4.0*g0);
  G4double A1 = std::max(0.0,(A0*g0 + A*(A-2*P-1)*0.25)/g1); 

  G4double E0 = U - A0;
  if (E0 <= 0.0) { return 0.0; }

  G4double E1 = std::max(0.0,theMaxKinEnergy - eKin - A1); 

  G4double Aj = A*(A+1)/(4.0*gj); 
  G4double Ej = std::max(0.0,efinal - Aj); 

  G4double rj = GetRj(P, aFragment.GetNumberOfCharged());
  G4double xs = CrossSection(eKin);
  
  G4double pA = fact*eKin*xs*rj 
    * CoalescenceFactor(theFragA) * FactorialFactor(N,P)
    * std::sqrt(2.0/(theReducedMass*efinal)) 
    * g4calc->powN(g1*E1/(g0*E0), N-A-1)
    * g4calc->powN(gj*Ej/(g0*E0), A-1)*gj*g1/(g0*g0*E0*theResA); 
   
  return pA;
}

G4double G4PreCompoundIon::GetBeta() const
{
  return -theCoulombBarrier;
}
