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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PreCompoundNeutron.cc,v 1.10 2001/11/22 17:38:09 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// by V. Lara

#include "G4PreCompoundNeutron.hh"


G4double G4PreCompoundNeutron::ProbabilityDistributionFunction(const G4double & eKin,
							       const G4Fragment & aFragment)
{
    if ( (aFragment.GetNumberOfParticles()-aFragment.GetNumberOfCharged()) < 1) 
	return 0.0;
    
    const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();
    // g = 0.595*a*A
    const G4double g = 0.595*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*GetRestA(); 
    G4double Alpha = 0.76+2.2/pow(GetRestA(),1.0/3.0); 
    G4double Beta = (2.12/pow(GetRestA(),2.0/3.0)-0.05)*MeV/Alpha;
	
    G4double Probability = 2.0/(pi*hbarc*hbarc*hbarc) * GetReducedMass() * Alpha * 
	r0 * r0 * pow(GetRestA(),2.0/3.0) *
	GetExcitonLevelDensityRatio()/(g*aFragment.GetExcitationEnergy()) *
	pow((1.0 - (eKin+GetBindingEnergy())/aFragment.GetExcitationEnergy()),
	    (aFragment.GetNumberOfExcitons()-2.0))*(eKin + Beta); 
	
    return Probability;
}

