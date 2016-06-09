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
// $Id: G4StatMFMacroTetraNucleon.cc,v 1.7 2008-10-24 22:57:07 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    const G4double BindingE = G4NucleiProperties::GetBindingEnergy(theA,2); //old value was 30.11*MeV
	
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
									
    return _Energy  = -G4NucleiProperties::GetBindingEnergy(theA,2) + 
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
