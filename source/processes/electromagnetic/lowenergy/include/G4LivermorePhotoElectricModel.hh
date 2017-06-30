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
// $Id: G4LivermorePhotoElectricModel.hh 104410 2017-05-30 07:17:09Z gcosmo $
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPhotoElectric developed by A.Forti and M.G.Pia
//
// 15 Mar 2010   L. Pandola, removed methods to set explicitely fluorescence cuts.
//               Main cuts from G4ProductionCutsTable are always used
// 30 May 2011   A Mantero & V Ivanchenko Migration to model design for deexcitation
// 22 Oct 2012   A & V Ivanchenko Migration data structure to G4PhysicsVector
//


#ifndef G4LivermorePhotoElectricModel_h
#define G4LivermorePhotoElectricModel_h 1

#include "G4VEmModel.hh"
#include "G4ElementData.hh"
#include <vector>

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;
class G4LPhysicsFreeVector;

class G4LivermorePhotoElectricModel : public G4VEmModel
{

public:

  G4LivermorePhotoElectricModel(const G4String& nam = "LivermorePhElectric");

  virtual ~G4LivermorePhotoElectricModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double energy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double energy, 
                                      G4double Z, 
                                      G4double A=0., 
                                      G4double cut=0.,
                                      G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);


  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

  inline void SetLimitNumberOfShells(G4int);

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:

  void ReadData(G4int Z, const char* path = nullptr);

  G4LivermorePhotoElectricModel & operator=(const G4LivermorePhotoElectricModel &right);
  G4LivermorePhotoElectricModel(const G4LivermorePhotoElectricModel&);

  G4ParticleDefinition*   theGamma;
  G4ParticleDefinition*   theElectron;

  G4int                   verboseLevel;
  G4int                   maxZ;
  G4int                   nShellLimit;
  G4bool                  fDeexcitationActive;
  G4bool                  isInitialised;

  static G4LPhysicsFreeVector*   fCrossSection[99];
  static G4LPhysicsFreeVector*   fCrossSectionLE[99];
  static std::vector<G4double>*  fParam[99];
  static G4int                   fNShells[99];
  static G4int                   fNShellsUsed[99];
  static G4ElementData*          fShellCrossSection;
  static G4Material*             fWater;
  static G4double                fWaterEnergyLimit;

  G4VAtomDeexcitation*    fAtomDeexcitation;

  G4double                fCurrSection;
  std::vector<G4double>   fSandiaCof;
};

inline 
void G4LivermorePhotoElectricModel::SetLimitNumberOfShells(G4int n)
{
  nShellLimit = n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
