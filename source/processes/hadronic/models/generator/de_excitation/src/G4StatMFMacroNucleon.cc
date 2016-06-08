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
// $Id: G4StatMFMacroNucleon.cc,v 1.7 2001/08/01 17:05:34 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroNucleon.hh"

// Operators

G4StatMFMacroNucleon & G4StatMFMacroNucleon::
operator=(const G4StatMFMacroNucleon & right)
{
    G4Exception("G4StatMFMacroNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroNucleon::operator==(const G4StatMFMacroNucleon & right) const
{
    G4Exception("G4StatMFMacroNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroNucleon::operator!=(const G4StatMFMacroNucleon & right) const
{
    G4Exception("G4StatMFMacroNucleon::operator!= meant to not be accessable");
    return true;
}

G4double G4StatMFMacroNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
						    const G4double nu, const G4double T)
{
    if (T <= 0.0) G4Exception("G4StatMFMacroNucleon::CalcMeanMultiplicity: Temperature less or equal 0");
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double degeneracy = 2.0;
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    _NeutronMeanMultiplicity = (degeneracy*FreeVol/lambda3)*exp(mu/T);
	
    _ProtonMeanMultiplicity = (degeneracy*FreeVol/lambda3)*
	exp((mu+nu-Coulomb)/T);

	

    return _MeanMultiplicity = _NeutronMeanMultiplicity + _ProtonMeanMultiplicity;
	
}


G4double G4StatMFMacroNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    return _Energy = Coulomb * theZARatio + (3./2.) * T;
							
}

G4double G4StatMFMacroNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double NeutronEntropy = 0.0;
    if (_NeutronMeanMultiplicity > 0.0)
	NeutronEntropy = _NeutronMeanMultiplicity*(5./2.+
						   log(2.0*G4double(theA)*FreeVol/(lambda3*_NeutronMeanMultiplicity)));
								
								
    G4double ProtonEntropy = 0.0;
    if (_ProtonMeanMultiplicity > 0.0)
	ProtonEntropy = _ProtonMeanMultiplicity*(5./2.+
						 log(2.0*G4double(theA)*FreeVol/(lambda3*_ProtonMeanMultiplicity)));
								
								
    return NeutronEntropy+ProtonEntropy;
}

