
#include "G4StatMFParameters.hh"


const G4double G4StatMFParameters::_Kappa = 1.0; // dimensionless

const G4double G4StatMFParameters::_KappaCoulomb = 2.0; // dimensionless

const G4double G4StatMFParameters::_Epsilon0 = 16.0*MeV;

// Bethe-Weizsacker coefficients
const G4double G4StatMFParameters::_E0 = 16.0*MeV;

const G4double G4StatMFParameters::_Beta0 = 18.0*MeV;

const G4double G4StatMFParameters::_Gamma0 = 25.0*MeV;

// Critical temperature (for liquid-gas phase transitions)
const G4double G4StatMFParameters::_CriticalTemp = 18.0*MeV;

// Nuclear radius
const G4double G4StatMFParameters::_r0 = 1.17*fermi;

G4double G4StatMFParameters::Beta(const G4double T)
{
	if (T > _CriticalTemp) return 0.0;
	else {
		G4double CriticalTempSqr = _CriticalTemp*_CriticalTemp;
		G4double TempSqr = T*T;
		G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
		return _Beta0*tmp*pow(tmp,1.0/4.0);
	}
}

G4double G4StatMFParameters::DBetaDT(const G4double T) 
{
	if (T > _CriticalTemp) return 0.0;
	else {
		G4double CriticalTempSqr = _CriticalTemp*_CriticalTemp;
		G4double TempSqr = T*T;
		G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
		return -5.0*_Beta0*pow(tmp,1.0/4.0)*(CriticalTempSqr*T)/
				((CriticalTempSqr+TempSqr)*(CriticalTempSqr+TempSqr));
	}
}


G4StatMFParameters G4StatMFParameters::theStatMFParameters;


G4StatMFParameters * G4StatMFParameters::GetAddress()
{ return &theStatMFParameters; }

