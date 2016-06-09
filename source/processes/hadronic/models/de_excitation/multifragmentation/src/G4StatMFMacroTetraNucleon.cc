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
// $Id: G4StatMFMacroTetraNucleon.cc,v 1.4 2005/06/04 13:27:49 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroTetraNucleon.hh"

// Copy constructor
G4StatMFMacroTetraNucleon::
G4StatMFMacroTetraNucleon(const G4StatMFMacroTetraNucleon & ) :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroTetraNucleon::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMacroTetraNucleon & G4StatMFMacroTetraNucleon::
operator=(const G4StatMFMacroTetraNucleon & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroTetraNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroTetraNucleon::operator==(const G4StatMFMacroTetraNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroTetraNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroTetraNucleon::operator!=(const G4StatMFMacroTetraNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroTetraNucleon::operator!= meant to not be accessable");
    return true;
}



G4double G4StatMFMacroTetraNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
							 const G4double nu, const G4double T)
{
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double degeneracy = 1;  // He4
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    const G4double BindingE = G4NucleiPropertiesTable::GetBindingEnergy(2,theA); //old value was 30.11*MeV
	
    G4double exponent = (BindingE + theA*(mu+nu*theZARatio+T*T/_InvLevelDensity) - 
			 Coulomb*theZARatio*theZARatio*std::pow(static_cast<G4double>(theA),5./3.))/T;
    if (exponent > 700.0) exponent = 700.0;
    
    _MeanMultiplicity = ( degeneracy*FreeVol* static_cast<G4double>(theA)* 
			  std::sqrt(static_cast<G4double>(theA))/lambda3)* 
	std::exp(exponent);
			 
    return _MeanMultiplicity;	
}


G4double G4StatMFMacroTetraNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    return _Energy  = -G4NucleiPropertiesTable::GetBindingEnergy(2,theA) + 
	Coulomb * theZARatio * theZARatio * std::pow(static_cast<G4double>(theA),5./3.) +
	(3./2.) * T +
	theA * T*T/_InvLevelDensity;
							
}



G4double G4StatMFMacroTetraNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double Entropy = 0.0;
    if (_MeanMultiplicity > 0.0)
	Entropy = _MeanMultiplicity*(5./2.+
				     std::log(8.0*FreeVol/(lambda3*_MeanMultiplicity)))+ // 8 = theA*std::sqrt(theA)
	    8.0*T/_InvLevelDensity;			
								
    return Entropy;
}
