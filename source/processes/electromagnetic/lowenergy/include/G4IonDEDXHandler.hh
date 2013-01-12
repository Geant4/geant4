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
// ===========================================================================
// GEANT4 class header file
//
// Class:                G4IonDEDXHandler
//
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 11. 03. 2009
//
// Modifications: 
//
//
// Class description:
//    Ion dE/dx table handler. 
//
// Comments:
//
// =========================================================================== 


#ifndef G4IONDEDXHANDLER_HH
#define G4IONDEDXHANDLER_HH

#include "globals.hh"
#include <vector>
#include <utility>
#include <list>
#include <map>

class G4ParticleDefinition;
class G4Material;
class G4PhysicsVector;
class G4VIonDEDXTable;
class G4VIonDEDXScalingAlgorithm;


// #########################################################################
// # Type definitions for a local cache
// #########################################################################

   typedef struct CacheValue{
      G4double energyScaling;         // Scaling factor for kinetic energy
      G4PhysicsVector* dedxVector;    // dE/dx vector for current projectile-
                                      // material combination
      G4double lowerEnergyEdge;       // Lower energy edge of dE/dx vector
      G4double upperEnergyEdge;       // Upper energy edge of dE/dx vector
      G4double density;               // Material density
   } G4CacheValue;


// #########################################################################
// # Class G4IonDEDXHandler: Handler class for stopping power tables
// #########################################################################

class G4IonDEDXHandler {

 public:
   G4IonDEDXHandler(G4VIonDEDXTable* tables,
                    G4VIonDEDXScalingAlgorithm* algorithm,
                    const G4String& name,
                    G4int maxCacheSize = 5,
                    G4bool splines = true);
   ~G4IonDEDXHandler();

   // Function checking the availability of stopping power values for a 
   // given ion-target combination (kinetic energy not considered) 
   G4bool IsApplicable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*);           // Target material             

   // Function returning the stopping power of a given material for a
   // projectile of specified energy
   G4double GetDEDX(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*,            // Target material 
              G4double);                    // Kinetic energy of projectile


   // Function for building stopping power vectors according to Bragg's
   // additivity rule
   G4bool BuildDEDXTable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*);           // Target material 

   // Function for building stopping power vectors according to Bragg's
   // additivity rule
   G4bool BuildDEDXTable(
              G4int atomicNumberIon,        // Atomic number of ion 
              const G4Material*);           // Target material 

   // Function printing stopping powers for a given ion-material combination
   // within a specified energy range 
   void PrintDEDXTable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material* ,           // Target material
              G4double,                     // Minimum energy per nucleon
              G4double,                     // Maximum energy per nucleon
              G4int,                        // Number of bins
              G4bool logScaleEnergy = true);// Logarithmic scaling of energy

   // Function returning the lower energy edge of stopping power tables
   G4double GetLowerEnergyEdge(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*);           // Target material 

   // Function returning the upper energy edge of stopping power tables 
   G4double GetUpperEnergyEdge(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*);           // Target material 

   // Function for clearing the cache
   void ClearCache();

   G4String GetName();

 private: 
   // The assignment operator and the copy constructor are hidden
   G4IonDEDXHandler& operator=(const G4IonDEDXHandler &r);
   G4IonDEDXHandler(const G4IonDEDXHandler&);

   // ######################################################################
   // # Stopping power table (table of stopping power vectors either built
   // # by G4VIonDEDXTable, or by the current class (using the Bragg 
   // # addivity rule)
   // ######################################################################

   // Class which creates dE/dx vectors 
   G4VIonDEDXTable* table;

   // Algorithm for scaling dE/dx values
   G4VIonDEDXScalingAlgorithm* algorithm;

   // Name associated with the dE/dx table
   G4String tableName;

   // Map of all dE/dx vectors
   typedef std::pair<G4int, const G4Material*> G4IonKey;
   typedef std::map<G4IonKey, G4PhysicsVector*> DEDXTable; 
   DEDXTable stoppingPowerTable;

   // Map of dE/dx vectors, built according to the Bragg additivity rule
   typedef std::map<G4IonKey, G4PhysicsVector*> DEDXTableBraggRule;
   DEDXTableBraggRule stoppingPowerTableBragg;

   // Flag indicating the usage of splines for dE/dx vectors built according
   // to Bragg rule
   G4bool useSplines;

   // ######################################################################
   // # "Most-recently-used" cache, to provide a faster access to physics 
   // # vectors
   // ######################################################################

   // A type definition of cache entry containing a key-value pair
   typedef std::pair<const G4ParticleDefinition*, const G4Material*> G4CacheKey;
   typedef struct CacheEntry {
      G4CacheKey key;
      G4CacheValue value;
   } G4CacheEntry;

   // A cache entry list, and a map of pointers to list iterators (for faster
   // searching)
   typedef std::list<G4CacheEntry> CacheEntryList;
   CacheEntryList cacheEntries;

   typedef std::map<G4CacheKey, void*> CacheIterPointerMap;
   CacheIterPointerMap cacheKeyPointers;  

   // Maximum number of cache entries
   G4int maxCacheEntries;

   // Function for updating the cache
   G4CacheValue UpdateCacheValue(
                const G4ParticleDefinition*,  // Projectile (ion) 
                const G4Material*);           // Target material

   // Function for retrieving cache values
   G4CacheValue GetCacheValue(
                const G4ParticleDefinition*,  // Projectile (ion) 
                const G4Material*);           // Target material
};

#endif  // G4IONDEDXHANDLER_HH
