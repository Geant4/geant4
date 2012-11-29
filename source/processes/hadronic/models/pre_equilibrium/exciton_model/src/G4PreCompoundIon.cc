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
// 10.02.2009 J. M. Quesada fixed bug in density level of light fragments  
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
  G4double r0 = theParameters->Getr0();
  fact = 0.75*CLHEP::millibarn/(CLHEP::pi*r0*r0*r0);
}

G4PreCompoundIon::~G4PreCompoundIon()
{}

G4double G4PreCompoundIon::
ProbabilityDistributionFunction(G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) { return 0.0; }
  G4double efinal = eKin + GetBindingEnergy();
  //G4cout << "Efinal= " << efinal << " Ekin= " << eKin << G4endl;
  if(efinal <= 0.0 ) { return 0.0; } 

  G4double U = aFragment.GetExcitationEnergy();
  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();
  G4int A = GetA();
  G4int N = P + H;

  G4double g0 = (6.0/pi2)*aFragment.GetA_asInt()*theParameters->GetLevelDensity();
  G4double g1 = (6.0/pi2)*GetRestA()*theParameters->GetLevelDensity();

  //JMQ 06/02/209  This is  THE BUG that was killing cluster emission
  // G4double gj = (6.0/pi2)*GetA() *
  //   G4PreCompoundParameters::GetAddress()->GetLevelDensity();

  G4double gj = g1;

  G4double A0 = G4double(P*P+H*H+P-3*H)/(4.0*g0);
  G4double A1 = std::max(0.0,(A0*g0 + A*(A-2*P-1)*0.25)/g1); 

  G4double E0 = U - A0;
  //G4cout << "E0= " << E0 << G4endl;
  if (E0 <= 0.0) { return 0.0; }

  G4double E1 = (std::max(0.0,GetMaximalKineticEnergy() - eKin - A1)); 

  G4double Aj = A*(A+1)/(4.0*gj); 
  G4double Ej = std::max(0.0,efinal - Aj); 

  G4double rj = GetRj(P, aFragment.GetNumberOfCharged());
  G4double xs = CrossSection(eKin);

  //G4cout << "rj= " << rj << " xs= " << xs << G4endl;

  // JMQ 10/02/09 reshaping of the formula (unnecessary std::pow elimitated)
  /*
  G4double r0 = theParameters->Getr0();
  G4double pA = (3.0/4.0) * std::sqrt(std::max(0.0, 2.0/(GetReducedMass()*
  		(eKin+GetBindingEnergy()))))/(pi * r0 * r0 *r0* GetRestA())* 
                eKin*CrossSection(eKin)*millibarn* 
                CoalescenceFactor(aFragment.GetA_asInt()) * FactorialFactor(N,P)* 
       GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged());

  G4double pB = std::pow((g1*E1)/(g0*E0),N-GetA()-1.0)*(g1/g0);
  G4double pC = std::pow((gj*Ej)/(g0*E0),GetA()-1.0)*(gj/g0)/E0; 
  pA *= pB * pC;
  */
  
  G4double pA = fact*eKin*xs*rj 
    * CoalescenceFactor(aFragment.GetA_asInt()) * FactorialFactor(N,P)
    * std::sqrt(2.0/(GetReducedMass()*efinal)) 
    * g4pow->powN(g1*E1/(g0*E0), N-A-1)
    * g4pow->powN(gj*Ej/(g0*E0), A-1)*gj*g1/(g0*g0*E0*GetRestA()); 
   
  return pA;
}
