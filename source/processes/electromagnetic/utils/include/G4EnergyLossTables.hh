// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossTables.hh,v 1.4 1999-04-13 09:01:11 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id:


#ifndef included_G4EnergyLossTables
#define included_G4EnergyLossTables

#include "globals.hh"
#include <rw/tvhdict.h>

#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ios.hh"

// A utility class, containing the energy loss tables
// for each particle
//
// Energy loss processes have to register their tables with this 
// class. The responsibility of creating and deleting the tables
// remains with the energy loss classes.
//
// P. Urban, 06/04/1998
// L. Urban, 27/05/1988 , modifications + new functions added
// L.Urban , 13/10/98 , revision
// L.Urban,  26/10/98 , revision, Interpolate removed 
// L.Urban , 08/02/99,  cache mechanism 
// L.Urban , 12/04/99 , bug fixed
// don't use the helper class.
// It can't be hidden for Rogue Wave uses it.

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
    G4Material *aMaterial);
  static G4double GetRange(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial);  
  static G4double GetLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial);  
  static G4double GetDeltaLabTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    G4Material *aMaterial);  
  static G4double GetProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial);  
  static G4double GetDeltaProperTime(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergyStart,
    G4double KineticEnergyEnd,
    G4Material *aMaterial);  

  static G4double GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial);
  static G4double GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial);  
  static G4double GetPreciseEnergyFromRange(
    const G4ParticleDefinition *aParticle,
    G4double range,
    G4Material *aMaterial); 
  
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
  static RWTValHashDictionary<K, G4EnergyLossTablesHelper> dict;
  
  static unsigned HashFun(const K& particle);

  static G4EnergyLossTablesHelper GetTables(const G4ParticleDefinition* p);

  static G4EnergyLossTablesHelper t ;
  static const G4ParticleDefinition* lastParticle ;
  static G4double Chargesquare ;
  static G4int oldIndex ;
  static G4double rmin,rmax,Thigh ;

};

#include "G4EnergyLossTables.icc"

#endif
