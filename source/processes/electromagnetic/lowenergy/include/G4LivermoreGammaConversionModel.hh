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
// Author: Sebastien Incerti
//         22 January 2012
//         on base of G4LivermoreGammaConversionModel (original version)
//         and G4LivermoreRayleighModel (MT version)

#ifndef G4LivermoreGammaConversionModel_h
#define G4LivermoreGammaConversionModel_h 1

#include "G4VEmModel.hh"
#include "G4Log.hh"

class G4ParticleChangeForGamma;
class G4LPhysicsFreeVector;
class G4PhysicsLogVector;

class G4LivermoreGammaConversionModel : public G4VEmModel
{

public:

  explicit G4LivermoreGammaConversionModel(
                      const G4ParticleDefinition* p = nullptr, 
		      const G4String& nam = "LivermoreConversion");

  virtual ~G4LivermoreGammaConversionModel();

  virtual void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			             G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0.0, 
                                      G4double cut=0.0,
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

  void ReadData(size_t Z, const char* path = nullptr);
  void InitialiseProbability(const G4ParticleDefinition*, G4int Z);

  inline G4double ScreenFunction1(G4double screenVariable);
  inline G4double ScreenFunction2(G4double screenVariable);

  G4LivermoreGammaConversionModel & operator=
  (const  G4LivermoreGammaConversionModel &right) = delete;
  G4LivermoreGammaConversionModel(const  G4LivermoreGammaConversionModel&) = delete;

  static G4double lowEnergyLimit;  
  static G4double tripletLowEnergy;
  static G4double tripletHighEnergy;
  
  static G4int verboseLevel;
  static G4int nbinsTriplet;
  static G4int maxZ;

  static G4LPhysicsFreeVector* data[100]; // 100 because Z range is 1-99
  static G4PhysicsLogVector*   probTriplet[100]; // 
  
  G4ParticleChangeForGamma* fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4LivermoreGammaConversionModel::ScreenFunction1(G4double screenVariable)
{
  // Compute the value of the screening function 3*phi1 - phi2
  return (screenVariable > 1.) 
    ? 42.24 - 8.368 * G4Log(screenVariable + 0.952)
    : 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4LivermoreGammaConversionModel::ScreenFunction2(G4double screenVariable)
{
  // Compute the value of the screening function 1.5*phi1 - 0.5*phi2
  return (screenVariable > 1.)
    ? 42.24 - 8.368 * G4Log(screenVariable + 0.952)
    : 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
