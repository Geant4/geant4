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
#include "G4UnitsTable.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1985Nuclear::G4hZiegler1985Nuclear():G4VhNuclearStoppingPower()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hZiegler1985Nuclear::~G4hZiegler1985Nuclear() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hZiegler1985Nuclear::NuclearStoppingPower(G4double kineticEnergy,
                                G4double z1, G4double z2, 
                                G4double m1, G4double m2) const
{  
  G4double energy = kineticEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (m1 + m2) * sqrt( pow(z1, .23) + pow(z2, .23) ) ;
  
  G4double er = 32.536 * m2 * energy / ( z1 * z2 * rm ) ;  // reduced energy
  
  if ( er <= 30 ) {
    ionloss = 0.5*log(1+1.1383*er)/
             (er+0.01312*pow(er,0.21226)+0.19593*sqrt(er)) ; 
    
  } else {
    ionloss = 0.5*log(er)/er ; 
  }
  
  // Stragling
  if(lossFlucFlag) {
    G4double sig = 4.0 * m1 * m2 / ((m1 + m2)*(m1 + m2)*
                  (4.0 + 0.197*pow(er,-1.6991)+6.584*pow(er,-1.0494))) ;


    ionloss *= G4RandGauss::shoot(1.0,sig) ;
  }

  ionloss *= 8.462 * z1 * z2 * m1 / rm ; // Return to [ev/(10^15 atoms/cm^2]

  if ( ionloss < 0.0) ionloss = 0.0 ;
  
  return ionloss;
}





