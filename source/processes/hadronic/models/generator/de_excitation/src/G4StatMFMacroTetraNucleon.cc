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
// $Id: G4StatMFMacroTetraNucleon.cc,v 1.6.2.1 2001/06/28 19:13:20 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroTetraNucleon.hh"

// Copy constructor
G4StatMFMacroTetraNucleon::
G4StatMFMacroTetraNucleon(const G4StatMFMacroTetraNucleon & right) :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    G4Exception("G4StatMFMacroTetraNucleon::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMacroTetraNucleon & G4StatMFMacroTetraNucleon::
operator=(const G4StatMFMacroTetraNucleon & right)
{
    G4Exception("G4StatMFMacroTetraNucleon::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMacroTetraNucleon::operator==(const G4StatMFMacroTetraNucleon & right) const
{
    G4Exception("G4StatMFMacroTetraNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroTetraNucleon::operator!=(const G4StatMFMacroTetraNucleon & right) const
{
    G4Exception("G4StatMFMacroTetraNucleon::operator!= meant to not be accessable");
    return true;
}



G4double G4StatMFMacroTetraNucleon::CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
							 const G4double nu, const G4double T)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
	
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
    const G4double degeneracy = 1;  // He4
	
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    const G4double BindingE = G4NucleiPropertiesTable::GetBindingEnergy(2,theA); //old value was 30.11*MeV
	
    _MeanMultiplicity = (degeneracy*FreeVol*G4double(theA)*sqrt(G4double(theA))/lambda3)* 
	exp((BindingE+ theA*(mu+nu*theZARatio+T*T/_InvLevelDensity) - 
	     Coulomb*theZARatio*theZARatio*pow(theA,5./3.))/T);
			 
    return _MeanMultiplicity;	
}


G4double G4StatMFMacroTetraNucleon::CalcEnergy(const G4double T)
{
    const G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
    return _Energy  = -G4NucleiPropertiesTable::GetBindingEnergy(2,theA) + 
	Coulomb * theZARatio * theZARatio * pow(G4double(theA),5./3.) +
	(3./2.) * T +
	theA * T*T/_InvLevelDensity;
							
}



G4double G4StatMFMacroTetraNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
    const G4double ThermalWaveLenght = 16.15*fermi/sqrt(T);
    const G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

    G4double Entropy = 0.0;
    if (_MeanMultiplicity > 0.0)
	Entropy = _MeanMultiplicity*(5./2.+
				     log(8.0*FreeVol/(lambda3*_MeanMultiplicity)))+ // 8 = theA*sqrt(theA)
	    8.0*T/_InvLevelDensity;			
								
    return Entropy;
}
