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
// Class:                G4IonParametrisedLossTable
// Helper classes:       G4IonLossTableHandle
//                       G4DummyScalingAlgorithm
// 
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 10. 11. 2008
//
// Modifications:
//
//
// Class descriptions:
//    G4IonParametrisedLossTable: Wrapper class for classes of the material
//    category, which maintain stopping power vectors for ions. This wrapper
//    class includes a cache for faster access of stopping powers and builds 
//    stopping power vectors for compounds using Bragg's additivity rule. 
//    G4IonLossTableHandle: A handle class for G4IonParametrisedLossTable
//    G4DummyScalingAlgorithm: A dummy algorithm to scale energies and
//    stopping powers (may be replaced with an algorithm, which relies on 
//    dynamic information of the particle: This may be required if stopping
//    power data is scaled using an effective charge approximation).
//
// Comments:
//
// =========================================================================== 


#ifndef G4IONPARAMETRISEDLOSSTABLE_HH
#define G4IONPARAMETRISEDLOSSTABLE_HH

#include "globals.hh"
#include "G4hIonEffChargeSquare.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>
#include <iomanip>
#include <map>
#include <list>
#include <utility>


// #########################################################################
// # Type definitions for a local cache
// # 
// #########################################################################

   typedef std::pair<const G4ParticleDefinition*, const G4Material*> G4CacheKey;

   typedef struct CacheValue{
      G4double energyScaling;         // Scaling factor for kinetic energy
      G4PhysicsVector* dedx;          // dEdx vector for current projectile-
                                      // material combination
      G4int dEdxIndex;                // dE/dx vector index
   } G4CacheValue;


// #########################################################################
// # Class G4IonLossTableHandle: Handle for table wrappers (to enable
// # vectors of table wrappers)   
// #########################################################################

class G4IonLossTableHandle {
  
 public:
   G4IonLossTableHandle() { }
   virtual ~G4IonLossTableHandle() { }

   virtual G4double GetLowerEnergyEdge(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*             // Target material 
	     ) = 0;
   virtual G4double GetUpperEnergyEdge(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*             // Target material 
	     ) = 0; 
   virtual G4bool IsApplicable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*             // Target material 
             ) = 0; 
   virtual G4bool IsApplicable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*,            // Target material 
              G4double                      // Kinetic energy of projectile
             ) = 0; 
   virtual G4double GetDEDX(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*,            // Target material 
              G4double                      // Kinetic energy of projectile
             ) = 0; 
   virtual void ClearCache() { }
};


// #########################################################################
// # Class G4DummyScalingAlgorithm:
// #########################################################################

class G4DummyScalingAlgorithm {

 public:
  // No scaling for kinetic energy
  inline G4double ScalingFactorEnergy(
            const G4ParticleDefinition*)    // Projectile (ion) 
                { return 1.0; }

  // No scaling for dE/dx
  inline G4double ScalingFactorDEDX(
            const G4ParticleDefinition*,    // Projectile (ion) 
            const G4Material*,              // Target material
            G4double)                       // Kinetic energy of projectile 
                { return 1.0; }
};


// #########################################################################
// # Class G4IonParametrisedLossTable: Wrapper for classes of the material
// # category, which maintain ion stopping power vectors
// #########################################################################

template <
    class LossTabulation,
    class ScalingAlgorithm = G4DummyScalingAlgorithm
>
class G4IonParametrisedLossTable : public G4IonLossTableHandle,
                                   public LossTabulation, 
                                   public ScalingAlgorithm {

 public:
   G4IonParametrisedLossTable(G4int maxSizeCache = 5,
                              G4bool splines = true);
   ~G4IonParametrisedLossTable();

   // Function checking the availability of stopping power values for a 
   // given ion-target combination and a given kinetic energy 
   G4bool IsApplicable(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*,            // Target material 
              G4double);                    // Kinetic energy of projectile

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
 
 private: 

   // Function for building stopping power vectors according to Bragg's
   // additivity rule
   G4PhysicsVector* BuildStoppingPowerVector(
              const G4ParticleDefinition*,  // Projectile (ion) 
              const G4Material*);           // Target material 

   // The assignment operator and the copy constructor are hidden
   G4IonParametrisedLossTable& operator=(const 
                                           G4IonParametrisedLossTable &right);
   G4IonParametrisedLossTable(const G4IonParametrisedLossTable&);

   // ######################################################################
   // # "Most-recently-used" cache, to provide a faster access to physics 
   // # vectors
   // ######################################################################

   // A type definition of cache entry containing a key-value pair
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

   // Function for clearing the cache
   void ClearCache();

   // ######################################################################
   // # Map of stopping power vectors computed according Bragg rule
   // #
   // ######################################################################

   typedef std::map<G4CacheKey, G4LPhysicsFreeVector*> StoppingPowerTable; 
   StoppingPowerTable s;

   // Flag indicating the usage of splines
   G4bool useSplines;

   // ######################################################################
   // # Energy boundaries of physics tables
   // ######################################################################

   G4double lowerEnergyEdge;
   G4double upperEnergyEdge;
};

#include "G4IonParametrisedLossTable.icc"

#endif  // G4IONPARAMETRISEDLOSSTABLE_HH
