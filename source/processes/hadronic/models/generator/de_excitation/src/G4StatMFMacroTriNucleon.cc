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
// $Id: G4StatMFMacroTriNucleon.cc,v 1.12 2003/05/30 13:23:26 hpw Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroTriNucleon.hh"

// Operators

G4StatMFMacroTriNucleon & G4StatMFMacroTriNucleon::
operator=(const G4StatMFMacroTriNucleon & )
{
    G4Exception("G4StatMFMacroTriNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroTriNucleon::operator==(const G4StatMFMacroTriNucleon & ) const
{
    G4Exception("G4StatMFMacroTriNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroTriNucleon::operator!=(const G4StatMFMacroTriNucleon & ) const
{
    G4Exception("G4StatMFMacroTriNucleon::operator!= meant to not be accessable");
    return true;
}



G4double G4StatMFMacroTriNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
						       const G4double nu, const G4double T)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double degeneracy = 2.0+2.0;  // H3 + He3
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    const G4double BindingE = G4NucleiPropertiesTable::GetBindingEnergy(1,theA); // old value was 9.224*MeV
//							+ G4NucleiProperties::GetBindingEnergy(2,theA);

    G4double exponent = (BindingE+ theA*(mu+nu*theZARatio) - 
			 Coulomb*theZARatio*theZARatio*pow(theA,5./3.))/T;
    if (exponent > 700.0) exponent = 700.0;

    _MeanMultiplicity = (degeneracy*FreeVol*static_cast<G4double>(theA)*
			 sqrt(static_cast<G4double>(theA))/lambda3)*
	exp(exponent);
			 
    return _MeanMultiplicity;
}


G4double G4StatMFMacroTriNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    return _Energy  = -G4NucleiPropertiesTable::GetBindingEnergy(1,theA) + 
	Coulomb * theZARatio * theZARatio * pow(static_cast<G4double>(theA),5./3.) +
	(3./2.) * T;
							
}


G4double G4StatMFMacroTriNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double Entropy = 0.0;
    if (_MeanMultiplicity > 0.0)
	Entropy = _MeanMultiplicity*(5./2.+
				     log(4.0*static_cast<G4double>(theA)*
					 sqrt(static_cast<G4double>(theA))*FreeVol/(lambda3*_MeanMultiplicity)));
								
								
    return Entropy;
}
