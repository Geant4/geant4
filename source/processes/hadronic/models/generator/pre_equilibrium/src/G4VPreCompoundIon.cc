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
// $Id: G4VPreCompoundIon.cc,v 1.8.2.1 2001/06/28 19:13:36 gunter Exp $
// GEANT4 tag $Name:  $
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
	

    // 	G4double R0J = 1.1;
    // 	G4double probA = GetCondensationProbability()*R0J*0.104/
    // 					(r0*pow(GetRestA(),1.0/3.0)*sqrt(GetA()*exEnergy));
    // 	G4double probB = GetExcitonLevelDensityRatio()*((eKin-GetCoulombBarrier())/exEnergy);
    // 	G4double ratio = (eKin+GetBindingEnergy())/exEnergy;
    // 	G4double exponent = GetRestA()-1.5;
    // 	if ( exponent>100. && ratio<1. ) return 0.;
    // 	G4double probC = pow( ratio, exponent );
    // 	G4double probD = pow( 1.0 - ratio,
    // 			aFragment.GetNumberOfExcitons()-GetA()-1.0 );
    // 	G4double prob = probA*probB*probC*probD;
	
	
    if (prob < 1.e-100) return 0.;
    else return prob;
}


G4double G4VPreCompoundIon::GetKineticEnergy(const G4Fragment & aFragment)
{
    G4double DJ = - GetCoulombBarrier();

    G4double T = aFragment.GetNumberOfParticles() + aFragment.GetNumberOfHoles() - GetA() - 1.0;
    G4double R2 = GetMaximalKineticEnergy();
    G4double R1 = R2 + GetCoulombBarrier();
	
    G4double E = 0.0;

    if (T <= -0.1) E = R1;
    else if (T <= 0.1) {
	G4double E1 = R1;
	G4double T3 = 0.0;
	do {
	    G4double PJ1 = GetA() - 1.5;
	    G4double AbsBindingE = abs(GetBindingEnergy());
	    if (GetBindingEnergy() <= 0.0 && AbsBindingE > GetCoulombBarrier()) {
		E = AbsBindingE + G4UniformRand()*aFragment.GetExcitationEnergy();
	    } else {
		E = GetCoulombBarrier() + G4UniformRand()*R2;
	    }
	    T3 = pow((E+GetBindingEnergy())/(E1+GetBindingEnergy()),PJ1)*((E+DJ)/(E1+DJ));
	} while (G4UniformRand() > T3);
    } else {
	G4double PJ1 = GetA() - 1.5;
	G4double ES = aFragment.GetExcitationEnergy()*(GetA()-0.5)+
	    (aFragment.GetExcitationEnergy()-R2)*(aFragment.GetNumberOfParticles()+
						  aFragment.GetNumberOfHoles()-2.5);
	G4double E1 = (ES + sqrt(ES*ES-(aFragment.GetExcitationEnergy()-R2)*(GetA()-1.5)*
				 (aFragment.GetNumberOfParticles()+aFragment.GetNumberOfHoles()-1.5)*
				 4.0*aFragment.GetExcitationEnergy()))/
	    ((aFragment.GetNumberOfParticles()+aFragment.GetNumberOfHoles()-1.5)*2.0)
	    - aFragment.GetExcitationEnergy() + R1;
	//
	G4double T3 = 0.0;
	do {
	    if (GetBindingEnergy() <= 0.0 && abs(GetBindingEnergy()) > GetCoulombBarrier()) {
		E = abs(GetBindingEnergy()) + G4UniformRand()*(aFragment.GetExcitationEnergy());
	    } else { 
		E = GetCoulombBarrier() + G4UniformRand()*R2;
	    }
	    T3 = (pow((E + GetBindingEnergy())/(E1 + GetBindingEnergy()),PJ1)*
		  ((E+DJ)/(E1+DJ))) * pow((R1-E)/(R1-E1),T);
	} while (G4UniformRand() > T3);
    }
    return E;
}
