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
// 
// -------------------------------------------------------------------
// $Id: G4PenelopeBremModel.hh,v 1.2 2004-07-27 08:38:30 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L.Pandola
//
// History:
// -----------
// 03 Feb 2003  L. Pandola       1st implementation
// 18 Mar 2003  L. Pandola       positrons added
// 23 May 2003  L. Pandola       ebins added as private member
// 25 May 2003  MGP              ebins removed from data members
// Class description:
// Penelope electromagnetic process, electron and positron BremModel
// --------------------------------------------------------------


#ifndef G4PENELOPEBremModel_HH
#define G4PENELOPEBremModel_HH 1

#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "globals.hh"
#include "G4PenelopeBremsstrahlungAngular.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;

class G4PenelopeBremModel : public G4VEmModel
{
  typedef std::vector<G4PenelopeBremsstrahlungAngular*> G4AngularData;
  //vector of pointers to the angular factors of the elements in each material

public:

  G4PenelopeBremModel(const G4ParticleDefinition* p=0,const G4String& nam = "PenelopeBrem");

  ~G4PenelopeBremModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition* p);

  G4double LowEnergyLimit(const G4ParticleDefinition* p);

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4Material*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy);

  G4double CrossSection(const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double,
                                      G4double) {return 0;};

  virtual G4double MaxSecondaryEnergy(
				const G4DynamicParticle* dynParticle);
protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				            G4double kineticEnergy);

private:

  // Hide copy constructor and assignment operator as private
  G4PenelopeBremModel(const G4PenelopeBremModel& );
  G4PenelopeBremModel& operator = (const G4PenelopeBremModel& right);

  void LoadAngularData();

  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Vector of pointers to the vectors containing tha angular data
  // one element = one material
  std::vector<G4AngularData*> materialAngularData;

  const G4ParticleDefinition* particle;
  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double minThreshold;
  G4int       verboseLevel;
  G4int       totBins;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
G4double G4PenelopeBremModel::MaxSecondaryEnergy(
				 const G4DynamicParticle* dynParticle)
{
  return dynParticle->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
G4double G4PenelopeBremModel::MaxSecondaryEnergy(
                                 const G4ParticleDefinition*,
    				       G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif
