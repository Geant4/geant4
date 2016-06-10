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
// $Id: G4StatMFMacroMultiNucleon.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 25.07.08 I.Pshenichnov (in collaboration with Alexander Botvina and Igor 
//          Mishustin (FIAS, Frankfurt, INR, Moscow and Kurchatov Institute, 
//          Moscow, pshenich@fias.uni-frankfurt.de) fixed computation of the 
//          symmetry energy 

#include "G4StatMFMacroMultiNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Default constructor
G4StatMFMacroMultiNucleon::
G4StatMFMacroMultiNucleon() :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::default_constructor meant to not be accessable");
}

// Copy constructor
G4StatMFMacroMultiNucleon::
G4StatMFMacroMultiNucleon(const G4StatMFMacroMultiNucleon & ) :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMacroMultiNucleon & G4StatMFMacroMultiNucleon::
operator=(const G4StatMFMacroMultiNucleon & ) 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroMultiNucleon::operator==(const G4StatMFMacroMultiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroMultiNucleon::operator!=(const G4StatMFMacroMultiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator!= meant to not be accessable");
    return true;
}



G4double G4StatMFMacroMultiNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu,
							 const G4double nu, const G4double T)
{
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double A23 = std::pow(static_cast<G4double>(theA),2./3.);
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
	
    G4double exponent = (mu + nu*theZARatio+ G4StatMFParameters::GetE0() + T*T/_InvLevelDensity 
			 - G4StatMFParameters::GetGamma0()*(1.0 - 2.0*theZARatio)*
			 (1.0 - 2.0*theZARatio))*theA
	- G4StatMFParameters::Beta(T)*A23 - Coulomb*theZARatio*theZARatio*A23*theA;
	
    exponent /= T;
	
    if (exponent > 30.0) exponent = 30.0;
	
    _MeanMultiplicity = std::max((FreeVol * static_cast<G4double>(theA) * 
				    std::sqrt(static_cast<G4double>(theA))/lambda3) *
				   std::exp(exponent),1.0e-30);
    return _MeanMultiplicity;	
}


G4double G4StatMFMacroMultiNucleon::CalcZARatio(const G4double nu)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    G4double den = 8.0*G4StatMFParameters::GetGamma0()+2.0*Coulomb*std::pow(static_cast<G4double>(theA),2./3.);
    G4double num = 4.0*G4StatMFParameters::GetGamma0()+nu;
	
    return theZARatio = num/den;
	

}



G4double G4StatMFMacroMultiNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
	
    const G4double A23 = std::pow(static_cast<G4double>(theA),2./3.);

    // Volume term 
    G4double EVol = static_cast<G4double>(theA) * (T*T/_InvLevelDensity - G4StatMFParameters::GetE0());
	
    // Symmetry term
    G4double ESym = static_cast<G4double>(theA) * G4StatMFParameters::GetGamma0() *(1. - 2.* theZARatio) * (1. - 2.* theZARatio);
	
    // Surface term
    G4double ESurf = A23*(G4StatMFParameters::Beta(T) - T*G4StatMFParameters::DBetaDT(T));
 
    // Coulomb term
    G4double ECoul = Coulomb*A23*static_cast<G4double>(theA)*theZARatio*theZARatio;
	
    // Translational term
    G4double ETrans = (3./2.)*T;
	
       
    return _Energy = EVol + ESurf + ECoul + ETrans + ESym;
}


G4double G4StatMFMacroMultiNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double Entropy = 0.0;
    if (_MeanMultiplicity > 0.0) {
	// Volume term
	G4double SV = 2.0*static_cast<G4double>(theA)*T/_InvLevelDensity;
		
	// Surface term
	G4double SS = -G4StatMFParameters::DBetaDT(T)*std::pow(static_cast<G4double>(theA),2./3.);
		
	// Translational term
	G4double ST = (5./2.)+std::log(FreeVol * std::sqrt(static_cast<G4double>(theA)) * 
				  static_cast<G4double>(theA)/(lambda3*_MeanMultiplicity));
		
		
	Entropy = _MeanMultiplicity*(SV + SS + ST);
    }
								
								
    return Entropy;
}
