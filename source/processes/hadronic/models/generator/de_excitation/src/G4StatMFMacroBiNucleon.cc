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
	
	_MeanMultiplicity = (degeneracy*FreeVol*G4double(theA)*sqrt(G4double(theA))/lambda3)*
			exp((BindingE+ theA*(mu+nu*theZARatio) - 
			 Coulomb*theZARatio*theZARatio*pow(theA,5./3.))/T);
			 
	return 	_MeanMultiplicity;
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
