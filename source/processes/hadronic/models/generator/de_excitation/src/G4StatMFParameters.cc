
#include "G4StatMFParameters.hh"


const G4double G4StatMFParameters::Kappa = 1.0;
const G4double G4StatMFParameters::KappaCoulomb = 2.0;
const G4double G4StatMFParameters::Epsilon0 = 16.0; // MeV
// Bethe-Weizsacker coefficients
const G4double G4StatMFParameters::E0 = 16.0; // MeV
const G4double G4StatMFParameters::Beta0 = 18.0; // MeV 
const G4double G4StatMFParameters::Gamma0 = 25.0; // MeV
// Critical temperature (for liquid-gas phase transitions)
const G4double G4StatMFParameters::CriticalTemp = 18.0; // MeV
// Nuclear radius
const G4double G4StatMFParameters::r0 = 1.17; // fm



G4StatMFParameters G4StatMFParameters::theStatMFParameters;


G4StatMFParameters * G4StatMFParameters::GetAddress()
{ return &theStatMFParameters; }

