//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4BetheHeitlerModel.hh,v 1.3 2005/05/12 11:06:42 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BetheHeitlerModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 19.04.2005
//
// Modifications:
//
// Class Description:
//
// Implementation of gamma convertion to e+e- in the field of a nucleus 
// 

// -------------------------------------------------------------------
//

#ifndef G4BetheHeitlerModel_h
#define G4BetheHeitlerModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForGamma;

class G4BetheHeitlerModel : public G4VEmModel
{

public:

  G4BetheHeitlerModel(const G4ParticleDefinition* p = 0, 
		      const G4String& nam = "Bethe-Heitler");
 
  virtual ~G4BetheHeitlerModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cut,
                                      G4double emax);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

private:

  G4double InitializeCrossSectionPerAtom(G4double energy, G4double Z);

  G4double ScreenFunction1(G4double ScreenVariable);

  G4double ScreenFunction2(G4double ScreenVariable);

  // hide assignment operator
  G4BetheHeitlerModel & operator=(const G4BetheHeitlerModel &right);
  G4BetheHeitlerModel(const  G4BetheHeitlerModel&);

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleDefinition*     thePositron;
  G4ParticleChangeForGamma* fParticleChange;
  G4PhysicsTable*           theCrossSectionTable; 

  G4double                  lowGammaEnergy;
  G4double                  highGammaEnergy;

  G4int                     nbins;
  size_t                    indexZ[120];
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4BetheHeitlerModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  G4bool b;
  size_t iz = indexZ[G4int(Z)];
  G4double x = (((*theCrossSectionTable)[iz]))->GetValue(energy, b);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4BetheHeitlerModel::ScreenFunction1(G4double ScreenVariable)

// compute the value of the screening function 3*PHI1 - PHI2

{
   G4double screenVal;

   if (ScreenVariable > 1.)
     screenVal = 42.24 - 8.368*std::log(ScreenVariable+0.952);
   else
     screenVal = 42.392 - ScreenVariable*(7.796 - 1.961*ScreenVariable);

   return screenVal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4BetheHeitlerModel::ScreenFunction2(G4double ScreenVariable)

// compute the value of the screening function 1.5*PHI1 - 0.5*PHI2

{
   G4double screenVal;

   if (ScreenVariable > 1.)
     screenVal = 42.24 - 8.368*std::log(ScreenVariable+0.952);
   else
     screenVal = 41.405 - ScreenVariable*(5.828 - 0.8945*ScreenVariable);

   return screenVal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
