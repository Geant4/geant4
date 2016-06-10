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
// $Id: G4DNABornExcitationModel1.hh 90057 2015-05-11 22:25:50Z matkara $
//

#ifndef G4DNABornExcitationModel1_h
#define G4DNABornExcitationModel1_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNACrossSectionDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4NistManager.hh"

class G4DNABornExcitationModel1: public G4VEmModel
{
public:
  G4DNABornExcitationModel1(const G4ParticleDefinition* p = 0,
                           const G4String& nam = "DNABornExcitationModel");

  virtual ~G4DNABornExcitationModel1();

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector& = *(new G4DataVector()));

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual G4double GetPartialCrossSection(const G4Material*,
                                          G4int level,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);

  inline void SelectStationary(G4bool input); 

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  G4bool statCode;

  // Water density table
  const std::vector<G4double>* fpMolWaterDensity;

  G4bool isInitialised;
  G4int verboseLevel;
  const G4ParticleDefinition* fParticleDefinition;

  double fLowEnergy;
  double fHighEnergy;
  G4String fTableFile;
  G4DNACrossSectionDataSet* fTableData;

  // Partial cross section
  G4int RandomSelect(G4double energy);
  
  G4DNAWaterExcitationStructure waterStructure;
   
  //
  G4DNABornExcitationModel1 & operator=(const  G4DNABornExcitationModel1 &right);
  G4DNABornExcitationModel1(const  G4DNABornExcitationModel1&);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNABornExcitationModel1::SelectStationary (G4bool input)
{ 
    statCode = input; 
}		 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
