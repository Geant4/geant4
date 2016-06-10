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
// Authors: G.Depaola & F.Longo

#ifndef G4LivermoreNuclearGammaConversionModel_h
#define G4LivermoreNuclearGammaConversionModel_h 1

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"

class G4LivermoreNuclearGammaConversionModel : public G4VEmModel
{

public:

  G4LivermoreNuclearGammaConversionModel(const G4ParticleDefinition* p = 0, 
		                  const G4String& nam = "LivermoreNuclearConversion");

  virtual ~G4LivermoreNuclearGammaConversionModel();

  virtual void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&);

  //MT
  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			             G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);
  //END MT

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0, 
                                      G4double cut=0,
                                      G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				       G4double tmin,
				       G4double maxEnergy);

  virtual G4double MinPrimaryEnergy(const G4Material*,
				    const G4ParticleDefinition*,
				    G4double);

private:

  void ReadData(size_t Z, const char* path = 0);

  G4double ScreenFunction1(G4double screenVariable);
  G4double ScreenFunction2(G4double screenVariable);

  G4LivermoreNuclearGammaConversionModel & operator=(const  G4LivermoreNuclearGammaConversionModel &right);
  G4LivermoreNuclearGammaConversionModel(const  G4LivermoreNuclearGammaConversionModel&);

  G4bool isInitialised;
  G4int verboseLevel;

  G4double lowEnergyLimit;  
  G4double smallEnergy;
  
  //MT
  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100]; // 100 because Z range is 1-99
                                          // in LivermoreRayleighModel, 101
					  //  because Z range is 1-100
  //END MT
  
  G4ParticleChangeForGamma* fParticleChange;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
