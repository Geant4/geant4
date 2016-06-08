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
// $Id: G4StatMFMacroBiNucleon.cc,v 1.9 2002/06/06 17:57:35 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroBiNucleon.hh"

// Operators

G4StatMFMacroBiNucleon & G4StatMFMacroBiNucleon::
operator=(const G4StatMFMacroBiNucleon & right)
{
    G4Exception("G4StatMFMacroBiNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroBiNucleon::operator==(const G4StatMFMacroBiNucleon & right) const
{
    G4Exception("G4StatMFMacroBiNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroBiNucleon::operator!=(const G4StatMFMacroBiNucleon & right) const
{
    G4Exception("G4StatMFMacroBiNucleon::operator!= meant to not be accessable");
    return true;
}


G4double G4StatMFMacroBiNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
						      const G4double nu, const G4double T)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
    
    const G4double degeneracy = 3.0;
    
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
    
    const G4double BindingE = G4NucleiPropertiesTable::GetBindingEnergy(1,theA); //old value was 2.796*MeV
    G4double exponent = (BindingE + theA*(mu+nu*theZARatio) - 
			 Coulomb*theZARatio*theZARatio*pow(theA,5./3.))/T;

    // To avoid numerical problems
    if (exponent < -700.0) exponent = -700.0;
    else if (exponent > 700.0) exponent = 700.0;

    _MeanMultiplicity = (degeneracy*FreeVol*G4double(theA)*sqrt(G4double(theA))/lambda3)*
	exp(exponent);
			 
    return _MeanMultiplicity;
}


G4double G4StatMFMacroBiNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    _Energy  = -G4NucleiPropertiesTable::GetBindingEnergy(1,theA) + 
	Coulomb * theZARatio * theZARatio * pow(G4double(theA),5./3.) +
	(3./2.) * T;
							
    return 	_Energy;				
}



G4double G4StatMFMacroBiNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double Entropy = 0.0;
    if (_MeanMultiplicity > 0.0)
	// Is this formula correct?
	Entropy = _MeanMultiplicity*(5./2.+
				     log(3.0*G4double(theA)*sqrt(G4double(theA))*FreeVol/
					 (lambda3*_MeanMultiplicity)));
								
								
    return Entropy;
}
