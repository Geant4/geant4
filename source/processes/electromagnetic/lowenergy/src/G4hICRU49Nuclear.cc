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
// File name:     G4hICRU49Nuclear
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
// ICRU Report N49, 1993. Moliere model.
// G.Moliere "Theorie der Streuung schneller geladener Teilchen I;
// Einzelstreuungam abbgeschirmten Coulomb-Feld" Z. f. Naturforsch, A2,
// 133 (1947).
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hICRU49Nuclear.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hICRU49Nuclear::G4hICRU49Nuclear():G4VhNuclearStoppingPower()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hICRU49Nuclear::~G4hICRU49Nuclear() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hICRU49Nuclear::NuclearStoppingPower(G4double kineticEnergy,
                                G4double z1, G4double z2, 
                                G4double m1, G4double m2) const
{  
  G4double energy = kineticEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (m1 + m2) * sqrt( pow(z1, .23) + pow(z2, .23) ) ;
  
  G4double er = 32.536 * m2 * energy / ( z1 * z2 * rm ) ;  // reduced energy
  
  static G4double a[104][2] = {
    1.0E+8,	5.831E-8,
    8.0E+7,	7.288E-8,
    6.0E+7,	9.719E-8,
    5.0E+7,	1.166E-7,
    4.0E+7,	1.457E-7,
    3.0E+7,	1.942E-7,
    2.0E+7,	2.916E-7,
    1.5E+7,	3.887E-7,

    1.0E+7,	5.833E-7,
    8.0E+6,	7.287E-7,
    6.0E+6,	9.712E-7,
    5.0E+6,	1.166E-6,
    4.0E+6,	1.457E-6,
    3.0E+6,	1.941E-6,
    2.0E+6,	2.911E-6,
    1.5E+6,	3.878E-6,

    1.0E+6,	5.810E-6,
    8.0E+5,	7.262E-6,
    6.0E+5,	9.663E-6,
    5.0E+5,	1.157E-5,
    4.0E+5,	1.442E-5,
    3.0E+5,	1.913E-5,
    2.0E+5,	2.845E-5,
    1.5E+5,	3.762E-5,

    1.0E+5,	5.554E-5,
    8.0E+4,	6.866E-5,
    6.0E+4,	9.020E-5,
    5.0E+4,	1.070E-4,
    4.0E+4,	1.319E-4,
    3.0E+4,	1.722E-4,
    2.0E+4,	2.499E-4,
    1.5E+4,	3.248E-4,

    1.0E+4,	4.688E-4,
    8.0E+3,	5.729E-4,
    6.0E+3,	7.411E-4,
    5.0E+3,	8.718E-4,
    4.0E+3,	1.063E-3,
    3.0E+3,	1.370E-3,
    2.0E+3,	1.955E-3,
    1.5E+3,	2.511E-3,

    1.0E+3,	3.563E-3,
    8.0E+2,	4.314E-3,
    6.0E+2,	5.511E-3,
    5.0E+2,	6.430E-3,
    4.0E+2,	7.756E-3,
    3.0E+2,	9.855E-3,
    2.0E+2,	1.375E-2,
    1.5E+2,	1.736E-2,

    1.0E+2,	2.395E-2,
    8.0E+1,	2.850E-2,
    6.0E+1,	3.552E-2,
    5.0E+1,	4.073E-2,
    4.0E+1,	4.802E-2,
    3.0E+1,	5.904E-2,
    1.5E+1,	9.426E-2,

    1.0E+1,	1.210E-1,
    8.0E+0,	1.377E-1,
    6.0E+0,	1.611E-1,
    5.0E+0,	1.768E-1,
    4.0E+0,	1.968E-1,
    3.0E+0,	2.235E-1,
    2.0E+0,	2.613E-1,
    1.5E+0,	2.871E-1,

    1.0E+0,	3.199E-1,
    8.0E-1,	3.354E-1,
    6.0E-1,	3.523E-1,
    5.0E-1,	3.609E-1,
    4.0E-1,	3.693E-1,
    3.0E-1,	3.766E-1,
    2.0E-1,	3.803E-1,
    1.5E-1,	3.788E-1,

    1.0E-1,	3.711E-1,
    8.0E-2,	3.644E-1,
    6.0E-2,	3.530E-1,
    5.0E-2,	3.444E-1,
    4.0E-2,	3.323E-1,
    3.0E-2,	3.144E-1,
    2.0E-2,	2.854E-1,
    1.5E-2,	2.629E-1,

    1.0E-2,	2.298E-1,
    8.0E-3,	2.115E-1,
    6.0E-3,	1.883E-1,
    5.0E-3,	1.741E-1,
    4.0E-3,	1.574E-1,
    3.0E-3,	1.372E-1,
    2.0E-3,	1.116E-1,
    1.5E-3,	9.559E-2,

    1.0E-3,	7.601E-2,
    8.0E-4,	6.668E-2,
    6.0E-4,	5.605E-2,
    5.0E-4,	5.008E-2,
    4.0E-4,	4.352E-2,
    3.0E-4,	3.617E-2,
    2.0E-4,	2.768E-2,
    1.5E-4,	2.279E-2,

    1.0E-4,	1.723E-2,
    8.0E-5,	1.473E-2,
    6.0E-5,	1.200E-2,
    5.0E-5,	1.052E-2,
    4.0E-5,	8.950E-3,
    3.0E-5,	7.246E-3,
    2.0E-5,	5.358E-3,
    1.5E-5,	4.313E-3,
    0.0,	3.166E-3
  };
  
  for (G4int i=1; i<104; i++)
    {
      if (er > a[i][0]) {
	ionloss = 
         (a[i][1]-a[i-1][1])*(er-a[i-1][0])/(a[i][0]-a[i-1][0])+a[i-1][1];   
	break;
      }
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





