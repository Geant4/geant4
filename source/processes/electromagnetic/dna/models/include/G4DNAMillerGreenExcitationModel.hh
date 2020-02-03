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

#ifndef G4DNAMillerGreenExcitationModel_h
#define G4DNAMillerGreenExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "Randomize.hh"
#include "G4NistManager.hh"

#include <deque>

class G4DNAMillerGreenExcitationModel : public G4VEmModel
{

public:

  G4DNAMillerGreenExcitationModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "DNAMillerGreenExcitationModel");

  virtual ~G4DNAMillerGreenExcitationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax);

  virtual G4double GetPartialCrossSection(const G4Material*,
                                          G4int /*level*/,
                                          const G4ParticleDefinition*,
                                          G4double /*kineticEnergy*/);

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

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  G4bool isInitialised;
  G4int verboseLevel;
  
  // Cross section

  // Partial cross section
  
  G4double PartialCrossSection(G4double energy,G4int level, const G4ParticleDefinition* particle);

  G4double Sum(G4double energy, const G4ParticleDefinition* particle);

  G4int RandomSelect(G4double energy, const G4ParticleDefinition* particle);

  G4int nLevels;

  G4DNAWaterExcitationStructure waterExcitation;
  
  G4double S_1s(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveCharge,
		G4double shellNumber);

  G4double S_2s(G4double t, 
		G4double energyTransferred,  
		G4double slaterEffectiveCharge, 
		G4double shellNumber);

  G4double S_2p(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveCharge, 
		G4double shellNumber);

  G4double R(G4double t, 
	     G4double energyTransferred, 
	     G4double slaterEffectiveCharge, 
	     G4double shellNumber);

  G4double kineticEnergyCorrection[4]; // 4 is the particle type index
  G4double slaterEffectiveCharge[3][4];
  G4double sCoefficient[3][4];

  //
   
  G4DNAMillerGreenExcitationModel & operator=(const  G4DNAMillerGreenExcitationModel &right);
  G4DNAMillerGreenExcitationModel(const  G4DNAMillerGreenExcitationModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNAMillerGreenExcitationModel::SelectStationary (G4bool input)
{ 
    statCode = input; 
}		 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
