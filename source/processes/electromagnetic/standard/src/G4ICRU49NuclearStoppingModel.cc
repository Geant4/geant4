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
// $Id: G4ICRU49NuclearStoppingModel.cc,v 1.2 2009-11-10 19:25:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4ICRU49NuclearStoppingModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 09.04.2008 from G4MuMscModel
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the model of ICRU'49 nuclear stopping

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ICRU49NuclearStoppingModel.hh"
#include "Randomize.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Step.hh"
#include "G4Pow.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU49NuclearStoppingModel::ad[] = {0.0};
G4double G4ICRU49NuclearStoppingModel::ed[] = {0.0};

using namespace std;

G4ICRU49NuclearStoppingModel::G4ICRU49NuclearStoppingModel(const G4String& nam) 
  : G4VEmModel(nam),lossFlucFlag(false)
{
  theZieglerFactor = eV*cm2*1.0e-15;
  g4pow = G4Pow::GetInstance();
  if(ad[0] == 0.0) InitialiseNuclearStopping();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU49NuclearStoppingModel::~G4ICRU49NuclearStoppingModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU49NuclearStoppingModel::Initialise(const G4ParticleDefinition*, 
					      const G4DataVector&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4ICRU49NuclearStoppingModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
						const G4MaterialCutsCouple*,
						const G4DynamicParticle*,
						G4double, G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4ICRU49NuclearStoppingModel::ComputeDEDXPerVolume(const G4Material* mat,
						   const G4ParticleDefinition* p,
						   G4double kinEnergy,
						   G4double)
{
  G4double nloss = 0.0;
  if(kinEnergy <= 0.0) return nloss; 

  // projectile
  G4double m1 = p->GetPDGMass();
  G4double z1 = std::fabs(p->GetPDGCharge()/eplus);

  if(kinEnergy*proton_mass_c2/m1 > z1*z1*100.*MeV) return nloss;

  // Projectile nucleus
  m1 /= amu_c2;

  //  loop for the elements in the material
  G4int numberOfElements = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* atomDensity  = mat->GetAtomicNumDensityVector();
 
  for (G4int iel=0; iel<numberOfElements; iel++) {
    const G4Element* element = (*theElementVector)[iel] ;
    G4double z2 = element->GetZ();
    G4double m2 = element->GetA()*mole/g ;
    nloss += (NuclearStoppingPower(kinEnergy, z1, z2, m1, m2))
           * atomDensity[iel] ;
  }
  nloss *= theZieglerFactor;
  return nloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ICRU49NuclearStoppingModel::NuclearStoppingPower(G4double kineticEnergy,
						   G4double z1, G4double z2,
						   G4double m1, G4double m2)
{
  G4double energy = kineticEnergy/keV ;  // energy in keV
  G4double nloss = 0.0;
  G4double z12 = z1*z2;
  G4int iz1 = G4int(z1);
  G4int iz2 = G4int(z2);
  
  G4double rm;
  if(iz1 > 1) rm = (m1 + m2) * ( g4pow->Z23(iz1) + g4pow->Z23(iz2) );
  else        rm = (m1 + m2) * g4pow->Z13(iz2);

  G4double er = 32.536 * m2 * energy / ( z12 * rm ) ;  // reduced energy

  if (er >= ed[0]) nloss = ad[0];
  else {
    // the table is inverse in energy
    for (G4int i=102; i>=0; --i)
    {
      if (er <= ed[i]) {
	nloss = (ad[i] - ad[i+1])*(er - ed[i+1])/(ed[i] - ed[i+1]) + ad[i+1];
	break;
      }
    }
  }

  // Stragling
  if(lossFlucFlag) {
    //    G4double sig = 4.0 * m1 * m2 / ((m1 + m2)*(m1 + m2)*
    //              (4.0 + 0.197*std::pow(er,-1.6991)+6.584*std::pow(er,-1.0494))) ;
    G4double sig = 4.0 * m1 * m2 / ((m1 + m2)*(m1 + m2)*
				    (4.0 + 0.197/(er*er) + 6.584/er));

    nloss *= G4RandGauss::shoot(1.0,sig) ;
    lossFlucFlag = false;
  }
   
  nloss *= 8.462 * z12 * m1 / rm ; // Return to [ev/(10^15 atoms/cm^2]

  if ( nloss < 0.0) nloss = 0.0 ;

  return nloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ICRU49NuclearStoppingModel::InitialiseNuclearStopping()
{
  const G4double nuca[104][2] = {
  { 1.0E+8, 5.831E-8},
  { 8.0E+7, 7.288E-8},
  { 6.0E+7, 9.719E-8},
  { 5.0E+7, 1.166E-7},
  { 4.0E+7, 1.457E-7},
  { 3.0E+7, 1.942E-7},
  { 2.0E+7, 2.916E-7},
  { 1.5E+7, 3.887E-7},

  { 1.0E+7, 5.833E-7},
  { 8.0E+6, 7.287E-7},
  { 6.0E+6, 9.712E-7},
  { 5.0E+6, 1.166E-6},
  { 4.0E+6, 1.457E-6},
  { 3.0E+6, 1.941E-6},
  { 2.0E+6, 2.911E-6},
  { 1.5E+6, 3.878E-6},

  { 1.0E+6, 5.810E-6},
  { 8.0E+5, 7.262E-6},
  { 6.0E+5, 9.663E-6},
  { 5.0E+5, 1.157E-5},
  { 4.0E+5, 1.442E-5},
  { 3.0E+5, 1.913E-5},
  { 2.0E+5, 2.845E-5},
  { 1.5E+5, 3.762E-5},

  { 1.0E+5, 5.554E-5},
  { 8.0E+4, 6.866E-5},
  { 6.0E+4, 9.020E-5},
  { 5.0E+4, 1.070E-4},
  { 4.0E+4, 1.319E-4},
  { 3.0E+4, 1.722E-4},
  { 2.0E+4, 2.499E-4},
  { 1.5E+4, 3.248E-4},

  { 1.0E+4, 4.688E-4},
  { 8.0E+3, 5.729E-4},
  { 6.0E+3, 7.411E-4},
  { 5.0E+3, 8.718E-4},
  { 4.0E+3, 1.063E-3},
  { 3.0E+3, 1.370E-3},
  { 2.0E+3, 1.955E-3},
  { 1.5E+3, 2.511E-3},

  { 1.0E+3, 3.563E-3},
  { 8.0E+2, 4.314E-3},
  { 6.0E+2, 5.511E-3},
  { 5.0E+2, 6.430E-3},
  { 4.0E+2, 7.756E-3},
  { 3.0E+2, 9.855E-3},
  { 2.0E+2, 1.375E-2},
  { 1.5E+2, 1.736E-2},

  { 1.0E+2, 2.395E-2},
  { 8.0E+1, 2.850E-2},
  { 6.0E+1, 3.552E-2},
  { 5.0E+1, 4.073E-2},
  { 4.0E+1, 4.802E-2},
  { 3.0E+1, 5.904E-2},
  { 1.5E+1, 9.426E-2},

  { 1.0E+1, 1.210E-1},
  { 8.0E+0, 1.377E-1},
  { 6.0E+0, 1.611E-1},
  { 5.0E+0, 1.768E-1},
  { 4.0E+0, 1.968E-1},
  { 3.0E+0, 2.235E-1},
  { 2.0E+0, 2.613E-1},
  { 1.5E+0, 2.871E-1},

  { 1.0E+0, 3.199E-1},
  { 8.0E-1, 3.354E-1},
  { 6.0E-1, 3.523E-1},
  { 5.0E-1, 3.609E-1},
  { 4.0E-1, 3.693E-1},
  { 3.0E-1, 3.766E-1},
  { 2.0E-1, 3.803E-1},
  { 1.5E-1, 3.788E-1},

  { 1.0E-1, 3.711E-1},
  { 8.0E-2, 3.644E-1},
  { 6.0E-2, 3.530E-1},
  { 5.0E-2, 3.444E-1},
  { 4.0E-2, 3.323E-1},
  { 3.0E-2, 3.144E-1},
  { 2.0E-2, 2.854E-1},
  { 1.5E-2, 2.629E-1},

  { 1.0E-2, 2.298E-1},
  { 8.0E-3, 2.115E-1},
  { 6.0E-3, 1.883E-1},
  { 5.0E-3, 1.741E-1},
  { 4.0E-3, 1.574E-1},
  { 3.0E-3, 1.372E-1},
  { 2.0E-3, 1.116E-1},
  { 1.5E-3, 9.559E-2},

  { 1.0E-3, 7.601E-2},
  { 8.0E-4, 6.668E-2},
  { 6.0E-4, 5.605E-2},
  { 5.0E-4, 5.008E-2},
  { 4.0E-4, 4.352E-2},
  { 3.0E-4, 3.617E-2},
  { 2.0E-4, 2.768E-2},
  { 1.5E-4, 2.279E-2},

  { 1.0E-4, 1.723E-2},
  { 8.0E-5, 1.473E-2},
  { 6.0E-5, 1.200E-2},
  { 5.0E-5, 1.052E-2},
  { 4.0E-5, 8.950E-3},
  { 3.0E-5, 7.246E-3},
  { 2.0E-5, 5.358E-3},
  { 1.5E-5, 4.313E-3},
  { 0.0, 3.166E-3}
  };

  for(G4int i=0; i<104; ++i) {
    ed[i] = nuca[i][0];
    ad[i] = nuca[i][1];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
