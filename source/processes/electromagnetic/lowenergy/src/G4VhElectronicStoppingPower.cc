// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
//
// File name:     G4VhElectronicStoppingPower
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Low energy hadrons/ions electronic stopping power parametrisation
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VhElectronicStoppingPower.hh" 
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VhElectronicStoppingPower::G4VhElectronicStoppingPower():
   theHeMassAMU(4.0026)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VhElectronicStoppingPower::~G4VhElectronicStoppingPower() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4VhElectronicStoppingPower::HeEffChargeSquare(
                                      const G4double z, 
                                      const G4double kineticEnergyHe) const
{
  // The aproximation of He effective charge from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  static G4double c[6] = {0.2865,  0.1266, -0.001429,
                          0.02402,-0.01135, 0.001475} ;

  G4double e = log( G4std::max( 1.0, kineticEnergyHe/(keV*GetHeMassAMU()))) ; 
  G4double x = c[0] ;
  G4double y = 1.0 ;
  for (G4int i=1; i<6; i++) {
    y *= e ;
    x += y * c[i] ;
  }

  G4double w = 7.6 -  e ;
  w = 1.0 + (0.007 + 0.00005*z) * exp( -w*w ) ;
  w = 4.0 * (1.0 - exp(-x)) * w * w ;

  return w;
}
