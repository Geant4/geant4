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
// Author: H. N. Tran (Ton Duc Thang University)
// p, H, He, He+ and He++ models are assumed identical
// NIMB 343, 132-137 (2015)
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#ifndef G4DNAIonElasticModel_h
#define G4DNAIonElasticModel_h 1

#include <map>
#include "G4DNACrossSectionDataSet.hh"
#include "G4VEmModel.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4DNAIonElasticModel : public G4VEmModel
{

public:

  G4DNAIonElasticModel (const G4ParticleDefinition* p = 0,
                        const G4String& nam ="DNAIonElasticModel");

  virtual
  ~G4DNAIonElasticModel ();

  virtual void
  Initialise (const G4ParticleDefinition* particuleDefinition,
              const G4DataVector&);

  virtual G4double
  CrossSectionPerVolume (const G4Material* material,
                         const G4ParticleDefinition* p, G4double ekin,
                         G4double emin, G4double emax);

  virtual void
  SampleSecondaries (std::vector<G4DynamicParticle*>*,
                     const G4MaterialCutsCouple*, const G4DynamicParticle*,
                     G4double tmin, G4double maxEnergy);

  void
  SetKillBelowThreshold (G4double threshold);

  G4double
  GetKillBelowThreshold ()
  {
    return killBelowEnergy;
  }

  inline void SelectStationary(G4bool input); 

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  G4bool statCode;

  // Water density table
  const std::vector<G4double>* fpMolWaterDensity;

  G4double killBelowEnergy;
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  G4bool isInitialised;
  G4int verboseLevel;

  G4double fParticle_Mass;

  // Cross section
  G4DNACrossSectionDataSet* fpTableData;

  // Final state

  G4double
  Theta (G4ParticleDefinition * aParticleDefinition, G4double k,
         G4double integrDiff);

  G4double
  LinLogInterpolate (G4double e1, G4double e2, G4double e, G4double xs1,
                     G4double xs2);

  G4double
  LogLogInterpolate (G4double e1, G4double e2, G4double e, G4double xs1,
                     G4double xs2);

  G4double
  QuadInterpolator (G4double e11, G4double e12, G4double e21, G4double e22,
                    G4double x11, G4double x12, G4double x21, G4double x22,
                    G4double t1, G4double t2, G4double t, G4double e);

  G4double
  LinLinInterpolate (G4double e1, G4double e2, G4double e, G4double xs1,
                     G4double xs2);

  typedef std::map<double, std::map<double, double> > TriDimensionMap;
  TriDimensionMap fDiffCrossSectionData;

  std::vector<double> eTdummyVec;

  typedef std::map<double, std::vector<double> > VecMap;
  VecMap eVecm;

  G4double
  RandomizeThetaCM (G4double k, G4ParticleDefinition * aParticleDefinition);

  //

  G4DNAIonElasticModel &
  operator= (const G4DNAIonElasticModel &right);
  G4DNAIonElasticModel (const G4DNAIonElasticModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNAIonElasticModel::SelectStationary (G4bool input)
{ 
    statCode = input; 
}		 

#endif

