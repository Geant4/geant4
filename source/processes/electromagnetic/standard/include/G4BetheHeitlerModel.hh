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
// $Id: G4BetheHeitlerModel.hh 104477 2017-06-01 07:39:33Z gcosmo $
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
// 02-02-06 Remove InitialiseCrossSectionPerAtom();
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
#include "G4Log.hh"

class G4ParticleChangeForGamma;
class G4Pow;

class G4BetheHeitlerModel : public G4VEmModel
{

public:

  explicit G4BetheHeitlerModel(const G4ParticleDefinition* p = 0, 
			       const G4String& nam = "BetheHeitler");
 
  virtual ~G4BetheHeitlerModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			       G4VEmModel* masterModel) override;

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0., 
                                      G4double cut=0.,
                                      G4double emax=DBL_MAX) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

private:

  G4double ScreenFunction1(G4double ScreenVariable);

  G4double ScreenFunction2(G4double ScreenVariable);

  // hide assignment operator
  G4BetheHeitlerModel & operator=(const G4BetheHeitlerModel &right) = delete;
  G4BetheHeitlerModel(const  G4BetheHeitlerModel&) = delete;

  G4Pow*                    g4calc;
  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleDefinition*     thePositron;
  G4ParticleChangeForGamma* fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4BetheHeitlerModel::ScreenFunction1(G4double ScreenVariable)

// compute the value of the screening function 3*PHI1 - PHI2
{
  return (ScreenVariable > 1.)
    ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    : 42.392 - ScreenVariable*(7.796 - 1.961*ScreenVariable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4BetheHeitlerModel::ScreenFunction2(G4double ScreenVariable)
// compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
{
  return (ScreenVariable > 1.)
    ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    : 41.405 - ScreenVariable*(5.828 - 0.8945*ScreenVariable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
