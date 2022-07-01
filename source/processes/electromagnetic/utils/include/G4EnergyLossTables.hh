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
//
//
// -----------------------------------------------------------------------------

#ifndef included_G4EnergyLossTables
#define included_G4EnergyLossTables

#include <map>
#include "globals.hh"

#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ios.hh"

//------------------------------------------------------------------------------
// A utility class, containing the energy loss tables
// for each particle
//
// Energy loss processes have to register their tables with this
// class. The responsibility of creating and deleting the tables
// remains with the energy loss classes.
// -----------------------------------------------------------------------------
//
// P. Urban, 06/04/1998
// L. Urban, 27/05/1988 , modifications + new functions added
// L.Urban , 13/10/98 , revision
// L.Urban,  26/10/98 , revision, Interpolate removed
// L.Urban , 08/02/99,  cache mechanism
// L.Urban , 12/04/99 , bug fixed
// don't use the helper class.
// It can't be hidden for Rogue Wave uses it.
// 10.11.99: moved from RWT hash dictionary to STL map, G.Barrand, M.Maire
// 26.10.01: all static functions movev from .icc to .cc file (mma)
// 15.01.03 Add interfaces required for "cut per region" (V.Ivanchenko)
// 12.03.03 Add warnings to obsolete interfaces (V.Ivanchenko)
// 12.04.03 move exception to new method (V.Ivanchenko)
//
// -----------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EnergyLossTablesHelper {

friend class G4EnergyLossTables;
  // the only instances are within the class G4EnergyLossTables

public:
  G4EnergyLossTablesHelper();

private:
  G4EnergyLossTablesHelper(const G4PhysicsTable* aDEDXTable,
			   const G4PhysicsTable* aRangeTable,
                           const G4PhysicsTable* anInverseRangeTable,
                           const G4PhysicsTable* aLabTimeTable,
                           const G4PhysicsTable* aProperTimeTable,
			   G4double aLowestKineticEnergy,
			   G4double aHighestKineticEnergy,
			   G4double aMassRatio,
                           G4int aNumberOfBins);
  // data to be stored in the dictionary
  const G4PhysicsTable* theDEDXTable;
  const G4PhysicsTable* theRangeTable;
  const G4PhysicsTable* theInverseRangeTable;
  const G4PhysicsTable* theLabTimeTable;
  const G4PhysicsTable* theProperTimeTable;
  G4double theLowestKineticEnergy;
  G4double theHighestKineticEnergy;
  G4double theMassRatio;
  G4int theNumberOfBins;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MaterialCutsCouple;

class G4EnergyLossTables {

public:

  // get the table for a given particle
  // (0 if the table was not found)
  static const G4PhysicsTable* GetDEDXTable(
    const G4ParticleDefinition* p);
  static const G4PhysicsTable* GetRangeTable(
    const G4ParticleDefinition* p);
  static const G4PhysicsTable* GetInverseRangeTable(
    const G4ParticleDefinition* p);
  static const G4PhysicsTable* GetLabTimeTable(
    const G4ParticleDefinition* p);
  static const G4PhysicsTable* GetProperTimeTable(
    const G4ParticleDefinition* p);

  // get the DEDX or the range for a given particle/energy/material
  static G4double GetDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetRange(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetDeltaLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    const G4Material *aMaterial);
  static G4double GetProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetDeltaProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    const G4Material *aMaterial);

  static G4double GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4Material *aMaterial);
  static G4double GetPreciseEnergyFromRange(
    const G4ParticleDefinition *aParticle,
    G4double range,
    const G4Material *aMaterial);

  // get the DEDX or the range for a given particle/energy/materialCutsCouple
  static G4double GetDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple,
    G4bool check = true);
  static G4double GetRange(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple,
    G4bool check = true);

  static G4double GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple);
  static G4double GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    const G4MaterialCutsCouple *couple);
  static G4double GetPreciseEnergyFromRange(
    const G4ParticleDefinition *aParticle,
    G4double range,
    const G4MaterialCutsCouple *couple,
    G4bool check = true);

  // to be called only by energy loss processes
  static void Register(
    const G4ParticleDefinition* p,
    const G4PhysicsTable* tDEDX,
    const G4PhysicsTable* tRange,
    const G4PhysicsTable* tInverseRange,
    const G4PhysicsTable* tLabTime,
    const G4PhysicsTable* tProperTime,
    G4double lowestKineticEnergy,
    G4double highestKineticEnergy,
    G4double massRatio,
    G4int NumberOfBins);

public:
  typedef const G4ParticleDefinition* K;

private:

  static void CPRWarning();
  static void ParticleHaveNoLoss(const G4ParticleDefinition* aParticle,
                                 const G4String&);

  /*
  typedef std::map<K,G4EnergyLossTablesHelper,std::less<K> > helper_map;
  static G4ThreadLocal helper_map *dict;

  static G4EnergyLossTablesHelper GetTables(const G4ParticleDefinition* p);

  static G4ThreadLocal G4EnergyLossTablesHelper *t;
  static G4ThreadLocal G4EnergyLossTablesHelper *null_loss;
  static G4ThreadLocal G4ParticleDefinition* lastParticle ;
  static G4ThreadLocal G4double QQPositron ;
  static G4ThreadLocal G4double Chargesquare ;
  static G4ThreadLocal G4int oldIndex ;
  static G4ThreadLocal G4double rmin,rmax,Thigh ;
  static G4ThreadLocal G4int  let_counter;
  static G4ThreadLocal G4int  let_max_num_warnings;
  static G4ThreadLocal G4bool first_loss;
  */
  using helper_map = std::map<K, G4EnergyLossTablesHelper, std::less<K> >;
  static helper_map *dict;

  static G4EnergyLossTablesHelper GetTables(const G4ParticleDefinition* p);

  static G4EnergyLossTablesHelper *t;
  static G4EnergyLossTablesHelper *null_loss;
  static G4ParticleDefinition* lastParticle ;
  static G4double QQPositron ;
  static G4double Chargesquare ;
  static G4int oldIndex ;
  static G4double rmin,rmax,Thigh ;
  static G4int  let_counter;
  static G4int  let_max_num_warnings;
  static G4bool first_loss;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
