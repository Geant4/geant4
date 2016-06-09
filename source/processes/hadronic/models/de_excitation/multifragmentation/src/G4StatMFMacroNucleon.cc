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
// $Id: G4StatMFMacroNucleon.cc,v 1.6 2008-07-25 11:20:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroNucleon.hh"

// Operators

G4StatMFMacroNucleon & G4StatMFMacroNucleon::
operator=(const G4StatMFMacroNucleon & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroNucleon::operator==(const G4StatMFMacroNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroNucleon::operator!=(const G4StatMFMacroNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroNucleon::operator!= meant to not be accessable");
    return true;
}

G4double G4StatMFMacroNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
						    const G4double nu, const G4double T)
{
    if (T <= 0.0) throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroNucleon::CalcMeanMultiplicity: Temperature less or equal 0");
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double degeneracy = 2.0;
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    G4double exponent_proton = (mu+nu-Coulomb)/T;
    G4double exponent_neutron = mu/T;

    if (exponent_neutron > 700.0) exponent_neutron = 700.0;
    if (exponent_proton > 700.0) exponent_proton = 700.0;

    _NeutronMeanMultiplicity = (degeneracy*FreeVol/lambda3)*std::exp(exponent_neutron);
	
    _ProtonMeanMultiplicity = (degeneracy*FreeVol/lambda3)*std::exp(exponent_proton);

	

    return _MeanMultiplicity = _NeutronMeanMultiplicity + _ProtonMeanMultiplicity;
	
}


G4double G4StatMFMacroNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    return _Energy = Coulomb * theZARatio * theZARatio + (3./2.) * T;
							
}

G4double G4StatMFMacroNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double NeutronEntropy = 0.0;
    if (_NeutronMeanMultiplicity > 0.0)
	NeutronEntropy = _NeutronMeanMultiplicity*(5./2.+
						   std::log(2.0*static_cast<G4double>(theA)*FreeVol/
						       (lambda3*_NeutronMeanMultiplicity)));
								
								
    G4double ProtonEntropy = 0.0;
    if (_ProtonMeanMultiplicity > 0.0)
	ProtonEntropy = _ProtonMeanMultiplicity*(5./2.+
						 std::log(2.0*static_cast<G4double>(theA)*FreeVol/
						     (lambda3*_ProtonMeanMultiplicity)));
								
								
    return NeutronEntropy+ProtonEntropy;
}

