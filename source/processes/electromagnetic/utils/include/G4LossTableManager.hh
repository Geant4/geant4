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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4LossTableManager
//
// Author:        Vladimir Ivanchenko on base of G4LossTables class
//                and Maria Grazia Pia ideas
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
//
// Class Description:
//
// A utility static class, responsable for the energy loss tables
// for each particle
//
// Energy loss processes have to register their tables with this
// class. The responsibility of creating and deleting the tables
// remains with the energy loss classes.

// -------------------------------------------------------------------
//

#ifndef G4LossTableManager_h
#define G4LossTableManager_h 1


#include "g4std/map"
#include "g4std/vector"
#include "globals.hh"
#include "G4LossTableBuilder.hh"
#include "G4VEnergyLossSTD.hh"

class G4PhysicsTable;
class G4MaterialCutsCouple;
class G4EnergyLossMessengerSTD;
class G4ParticleDefinition;

class G4LossTableManager
{

public:

  static G4LossTableManager* Instance();

  ~G4LossTableManager();

  void Clear();

  // get the DEDX or the range for a given particle/energy/material
  G4double GetDEDX(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4MaterialCutsCouple *couple);

  G4double GetRange(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4MaterialCutsCouple *couple);

  G4double GetEnergy(
    const G4ParticleDefinition *aParticle,
    G4double range,
    const G4MaterialCutsCouple *couple);

  // to be called only by energy loss processes
  void Register(G4VEnergyLossSTD* p);

  void BuildDEDXTable(const G4ParticleDefinition* aParticle);

  void RetrieveDEDXTable(const G4ParticleDefinition* aParticle,
                               G4VEnergyLossSTD*);

  void SetLossFluctuations(G4bool val);

  void SetSubCutoff(G4bool val);

  void SetIntegral(G4bool val);

  void SetRandomStep(G4bool val);

  void SetMinSubRange(G4double val);

  void SetMinEnergy(G4double val);

  void SetMaxEnergy(G4double val);

  void SetStepLimits(G4double v1, G4double v2);

  G4bool StorePhysicsTable(G4VEnergyLossSTD* el,
		     const G4String& dedx_file,
		     const G4String& range_file,
		     const G4String& inv_range_file,
			   G4bool ascii);

private:

  G4LossTableManager();

  void Initialise();

  G4VEnergyLossSTD* BuildTables(const G4ParticleDefinition* aParticle);

  //  G4LossTableManager(const G4LossTableManager&);
  //  const G4LossTableManager& operator=(const G4LossTableManager& right);

private:

  static G4LossTableManager* theInstance;

  typedef const G4ParticleDefinition* PD;
  G4std::map<PD,G4VEnergyLossSTD*,G4std::less<PD> > loss_map;

  G4std::vector<G4VEnergyLossSTD*> loss_vector;
  G4std::vector<PD> part_vector;
  G4std::vector<PD> base_part_vector;
  G4std::vector<G4bool> tables_are_built;
  G4std::vector<G4PhysicsTable*> dedx_vector;
  G4std::vector<G4PhysicsTable*> range_vector;
  G4std::vector<G4PhysicsTable*> inv_range_vector;

  // cash
  G4VEnergyLossSTD*    currentLoss;
  PD                   currentParticle;
  PD                   theElectron;

  G4int n_loss;

  G4bool all_tables_are_built;
  G4bool all_particles_are_mapped;
  G4bool electron_table_are_built;
  G4bool lossFluctuationFlag;
  G4bool subCutoffFlag;
  G4bool rndmStepFlag;
  G4bool integral;

  G4double minSubRange;
  G4double maxRangeVariation;
  G4double maxFinalStep;
  G4double minKinEnergy;
  G4double maxKinEnergy;

  G4VEnergyLossSTD*         eIonisation;
  G4LossTableBuilder*       tableBuilder;
  G4EnergyLossMessengerSTD* theMessenger;
  G4int verbose;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline G4double G4LossTableManager::GetDEDX(
          const G4ParticleDefinition *aParticle,
                G4double kineticEnergy,
          const G4MaterialCutsCouple *couple)
{
  if(aParticle != currentParticle) {
    G4std::map<PD, G4VEnergyLossSTD*,G4std::less<PD> >::const_iterator pos;
    if((pos = loss_map.find(aParticle)) == loss_map.end()) return 0.0;
    currentParticle = aParticle;
    currentLoss = (*pos).second;
  }
  return currentLoss->GetDEDX(kineticEnergy, couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline G4double G4LossTableManager::GetRange(
          const G4ParticleDefinition *aParticle,
                G4double kineticEnergy,
          const G4MaterialCutsCouple *couple)
{
  if(aParticle != currentParticle) {
    G4std::map<PD, G4VEnergyLossSTD*, G4std::less<PD> >::const_iterator pos;
    if((pos = loss_map.find(aParticle)) == loss_map.end()) return 0.0;
    currentParticle = aParticle;
    currentLoss = (*pos).second;
  }
  return currentLoss->GetRange(kineticEnergy, couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4LossTableManager::GetEnergy(
          const G4ParticleDefinition *aParticle,
                G4double range,
          const G4MaterialCutsCouple *couple)
{
  if(aParticle != currentParticle) {
    G4std::map<PD,G4VEnergyLossSTD*,G4std::less<PD> >::const_iterator pos;
    if((pos = loss_map.find(aParticle)) == loss_map.end()) return 0.0;
    currentParticle = aParticle;
    currentLoss = (*pos).second;
  }
  return currentLoss->GetKineticEnergy(range, couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

