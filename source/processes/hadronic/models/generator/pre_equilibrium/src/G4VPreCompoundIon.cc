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
// $Id: G4VPreCompoundIon.cc,v 1.12 2002/01/15 12:57:38 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// by V. Lara

#include "G4VPreCompoundIon.hh"  


G4double G4VPreCompoundIon::ProbabilityDistributionFunction(const G4double & eKin, 
							    const G4Fragment & aFragment)
{
    G4int pplus = aFragment.GetNumberOfCharged();
    G4int pneut = aFragment.GetNumberOfParticles()-pplus;
    if (pneut < (GetA()-GetZ()) || pplus < GetZ()) return 0.0;
    
    const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();
    G4double exEnergy = aFragment.GetExcitationEnergy();
    G4double B = GetBindingEnergy();
    
    G4double Z = aFragment.GetZ();
    G4double C = GetCCoef(Z);
    
    G4double probA = (3.0/4.0)*sqrt(2.0/GetReducedMass())*(1.0+C)*GetExcitonLevelDensityRatio()*
	GetCondensationProbability()*(eKin - GetCoulombBarrier())/
	(r0*pow(GetRestA(),1.0/3.0)*exEnergy*sqrt(eKin+B));
    
    G4double base = 1.0 + B/exEnergy;
    G4double exponent = GetA() - 1.0;
    if (exponent > 100.0 && base < 1.0) return 0.0;
    G4double probB = pow(base,exponent);
    
    base = 1.0 - ((eKin+B)/exEnergy);
    exponent = aFragment.GetNumberOfExcitons() - 1.0 - GetA();
    if (exponent > 100.0 && base < 1.0) return 0.0;
    G4double probC = pow(base,exponent);
    
    G4double prob = probA * probB * probC;
        
    if (prob < 1.e-100) return 0.;
    else return prob;
}


