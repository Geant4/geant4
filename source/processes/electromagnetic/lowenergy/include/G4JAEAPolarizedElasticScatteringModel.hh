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
/*/
Authors:
M. Omer and R. Hajima  on 15 November 2019
contact:
omer.mohamed@jaea.go.jp and hajima.ryoichi@qst.go.jp
Publication Information:
1- M. Omer, R. Hajima, Validating polarization effects in gamma-rays elastic scattering by Monte
Carlo simulation, New J. Phys., vol. 21, 2019, pp. 113006 (1-10),
https://doi.org/10.1088/1367-2630/ab4d8a
*/

#ifndef G4JAEAPolarizedElasticScatteringModel_h
#define G4JAEAPolarizedElasticScatteringModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DataVector.hh"

class G4JAEAPolarizedElasticScatteringModel : public G4VEmModel
{
public:
  explicit G4JAEAPolarizedElasticScatteringModel();
  virtual ~G4JAEAPolarizedElasticScatteringModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*,
		       G4VEmModel* masterModel) override;

  void InitialiseForElement(const G4ParticleDefinition*, G4int Z) override;

  G4double ComputeCrossSectionPerAtom(
				      const G4ParticleDefinition*,
                                      G4double kinEnergy,
                                      G4double Z,
                                      G4double A=0,
                                      G4double cut=0,
                                      G4double emax=DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  inline void SetLowEnergyThreshold(G4double val){lowEnergyLimit = val;};
  inline void SetPolarizationSensitvity(G4bool,G4bool,G4bool);
  
  void SetDebugVerbosity(G4int val){verboseLevel = val;};

  G4JAEAPolarizedElasticScatteringModel & operator=(const G4JAEAPolarizedElasticScatteringModel &right) = delete;
  G4JAEAPolarizedElasticScatteringModel(const G4JAEAPolarizedElasticScatteringModel&) = delete;

private:
  void ReadData(std::size_t Z, const char* path = 0);
  G4double GeneratePolarizedPhi(G4double Sigma_para,G4double Sigma_perp, G4double initial_Pol_Plane);

  static const G4int maxZ = 99;
  static G4PhysicsFreeVector* dataCS[maxZ+1];
  static G4DataVector* Polarized_ES_Data[maxZ+1];
  G4double distribution[181];
  G4double cdistribution[181];
    
  G4ParticleChangeForGamma* fParticleChange;
  G4double lowEnergyLimit;

  G4int verboseLevel;
  
  G4bool fLinearPolarizationSensitvity1;
  G4bool fLinearPolarizationSensitvity2;
  G4bool fCircularPolarizationSensitvity;
  G4bool isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4JAEAPolarizedElasticScatteringModel::SetPolarizationSensitvity(G4bool linear1,
									     G4bool linear2, 
									     G4bool circular)
{
  fLinearPolarizationSensitvity1=linear1;
  fLinearPolarizationSensitvity2=linear2;
  fCircularPolarizationSensitvity=circular;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
