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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
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
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Exp.hh"

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

  static const G4double c[6] = {0.2865,  0.1266, -0.001429,
				0.02402,-0.01135, 0.001475} ;

  G4double e = std::log( std::max( 1.0, kineticEnergyHe/(keV*GetHeMassAMU()))) ;
  G4double x = c[0] ;
  G4double y = 1.0 ;
  for (G4int i=1; i<6; i++) {
    y *= e ;
    x += y * c[i] ;
  }

  G4double w = 7.6 -  e ;
  w = 1.0 + (0.007 + 0.00005*z) * G4Exp( -w*w ) ;
  w = 4.0 * (1.0 - G4Exp(-x)) * w * w ;
  return w;
}
