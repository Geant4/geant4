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
// File name:     G4hZiegler1985Nuclear
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
// Nuclear stopping power parametrised according to
// J.F.Ziegler, J.P. Biersack, U. Littmark
// The Stopping and Range of Ions in Matter,
// Vol.1, Pergamon Press, 1985
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hZiegler1985Nuclear.hh"

#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1985Nuclear::G4hZiegler1985Nuclear():G4VhNuclearStoppingPower()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1985Nuclear::~G4hZiegler1985Nuclear() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hZiegler1985Nuclear::NuclearStoppingPower(G4double kineticEnergy,
                                G4double z1, G4double z2, 
                                G4double m1, G4double m2Local) const
{  
  G4double energy = kineticEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (m1 + m2Local) * ( std::pow(z1, .23) + std::pow(z2, .23))  ;

  G4double er = 32.536 * m2Local * energy / ( z1 * z2 * rm ) ;  // reduced energy

  if ( er <= 30 ) {
    ionloss = 0.5*std::log(1+1.1383*er)/
             (er+0.01312*std::pow(er,0.21226)+0.19593*std::sqrt(er)) ;

  } else {
    ionloss = 0.5*std::log(er)/er ;
  }
  
  // Stragling
  if(lossFlucFlag) {
    G4double sig = 4.0 * m1 * m2Local / ((m1 + m2Local)*(m1 + m2Local)*
                  (4.0 + 0.197*std::pow(er,-1.6991)+6.584*std::pow(er,-1.0494))) ;


    ionloss *= G4RandGauss::shoot(1.0,sig) ;
  }

  ionloss *= 8.462 * z1 * z2 * m1 / rm ; // Return to [ev/(10^15 atoms/cm^2]

  if ( ionloss < 0.0) ionloss = 0.0 ;
  
  return ionloss;
}





