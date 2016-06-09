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
//
// $Id: G4PreCompoundNucleon.cc,v 1.10.2.1 2009/03/03 13:17:04 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-03 $
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
//

#include "G4PreCompoundNucleon.hh"
#include "G4PreCompoundParameters.hh"

G4bool G4PreCompoundNucleon::IsItPossible(const G4Fragment& aFragment) 
{
  G4int pplus = aFragment.GetNumberOfCharged();   
  G4int pneut = aFragment.GetNumberOfParticles()-pplus;
  return (pneut >= (GetA()-GetZ()) && pplus >= GetZ());
}

G4double G4PreCompoundNucleon::
ProbabilityDistributionFunction(const G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) return 0.0;

  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  
  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
 
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();

  G4double A0 = ((P*P+H*H+P-H)/4.0 - H/2.0)/g0;

  G4double A1 = (A0 - P/2.0)/g1;
  
  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  G4double E1 = std::max(0.0,U - eKin - GetBindingEnergy() - A1);
  if (E1 == 0.0) return 0.0;


  G4double Probability = 1.0/pi2*2.0/(hbarc*hbarc*hbarc) * GetReducedMass()
    * GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())  
    * eKin*CrossSection(eKin)*millibarn * P*(N-1.0) * std::pow(g1*E1/(g0*E0),N-2.0)/E0 * g1/(g0*g0);

  if (GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())<0.0 
      || CrossSection(eKin) <0) {  
    std::ostringstream errOs;
    G4cout<<"WARNING:  NEGATIVE VALUES "<<G4endl;     
    errOs << "Rj=" << GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())
	  <<G4endl;
    errOs <<"  xsec("<<eKin<<" MeV) ="<<CrossSection(eKin)<<G4endl;
    errOs <<"  A="<<GetA()<<"  Z="<<GetZ()<<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }

  return Probability;
}
















