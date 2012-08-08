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
// File name:     G4hZiegler1977Nuclear
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
// J.F.Ziegler, Helium Stopping Powers and
// Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hZiegler1977Nuclear.hh"

#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1977Nuclear::G4hZiegler1977Nuclear():G4VhNuclearStoppingPower()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1977Nuclear::~G4hZiegler1977Nuclear() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hZiegler1977Nuclear::NuclearStoppingPower(G4double kineticEnergy,
                                G4double z1, G4double z2, 
                                G4double m1, G4double m2Local) const
{  
  G4double energy = kineticEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (m1 + m2Local) * std::sqrt( std::pow(z1, 0.667) + std::pow(z2, 0.667) ) ;
  
  G4double er = 32.53 * m2Local * energy / ( z1 * z2 * rm ) ;  // reduced energy
  
  if ( er < 0.01 ) {
    ionloss = std::sqrt(er) * 1.593 ; 
    
  } else if ( er < 10.0 ) {
    ionloss = 1.7 * std::sqrt(er) * std::log(er + std::exp(1.0)) / 
      (1.0 + 6.8 * er + 3.4 * std::pow(er, 1.5)) ; 
    
  } else {
    ionloss = std::log(0.47 * er) * 0.5 / er  ;
  }

  // Stragling
  if(lossFlucFlag) {
    G4double sig = 4.0 * m1 * m2Local * std::sqrt( (std::pow(z1, 0.23) + std::pow(z2, 0.23)) / 
                                         (std::pow(z1, 0.667) + std::pow(z2, 0.667)) ) 
                 / ((m1 +m2Local)*(m1 + m2Local)*
                    (4.0 + 0.197*std::pow(er,-1.6991)+6.584*std::pow(er,-1.0494))) ;

    ionloss *= G4RandGauss::shoot(1.0,sig) ;
  }
  
  ionloss *= 8.462 * z1 * z2 * m1 / rm ; // Return to [ev/(10^15 atoms/cm^2]

  if ( ionloss < 0.0) ionloss = 0.0 ;
  
  return ionloss;
}



