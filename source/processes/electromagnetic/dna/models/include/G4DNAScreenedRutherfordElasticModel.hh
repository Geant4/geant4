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

#ifndef G4DNAScreenedRutherfordElasticModel_h
#define G4DNAScreenedRutherfordElasticModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4DNAScreenedRutherfordElasticModel : public G4VEmModel
{
public:
  G4DNAScreenedRutherfordElasticModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "DNAScreenedRutherfordElasticModel");

  virtual ~G4DNAScreenedRutherfordElasticModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);

  inline void SetKillBelowThreshold (G4double threshold);

  inline void SelectFasterComputation(G4bool input); 

protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  G4double lowEnergyLimit;
  G4double intermediateEnergyLimit;
  G4double highEnergyLimit;
  
  // -- Brenner & Zaider
  std::vector<G4double> betaCoeff;
  std::vector<G4double> deltaCoeff;
  std::vector<G4double> gamma035_10Coeff;
  std::vector<G4double> gamma10_100Coeff;
  std::vector<G4double> gamma100_200Coeff;
  
  // -- Water density table
  const std::vector<G4double>* fpWaterDensity;
  
  G4int verboseLevel;
  
  G4bool isInitialised;
  G4bool fasterCode;
  
  // -- Cross section
  G4double RutherfordCrossSection(G4double energy, G4double z);
  G4double ScreeningFactor(G4double energy, G4double z);
  
  // -- Final state according to Brenner & Zaider
  G4double BrennerZaiderRandomizeCosTheta(G4double k);
  G4double CalculatePolynomial(G4double k, std::vector<G4double>& vec);

   
  // -- Final state according to Screened Rutherford
  G4double ScreenedRutherfordRandomizeCosTheta(G4double k, G4double z);

  //
  G4DNAScreenedRutherfordElasticModel& operator=
    (const G4DNAScreenedRutherfordElasticModel &right);
  G4DNAScreenedRutherfordElasticModel(const G4DNAScreenedRutherfordElasticModel&);
};

inline void
G4DNAScreenedRutherfordElasticModel::SetKillBelowThreshold (G4double)
{
    G4ExceptionDescription errMsg;
    errMsg << "The method G4DNAScreenedRutherfordElasticModel::"
              "SetKillBelowThreshold is deprecated";
    
    G4Exception("G4DNAScreenedRutherfordElasticModel::SetKillBelowThreshold",
                "deprecated",
                JustWarning,
                errMsg);
}		 

inline void
G4DNAScreenedRutherfordElasticModel::SelectFasterComputation(G4bool input)
{ 
    fasterCode = input; 
}		 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
