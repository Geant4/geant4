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
// GEANT4 class source file
//
// Class:                G4IonDEDXHandler
//
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 11. 03. 2009
//
// Modifications: 12. 11 .2009 - Function BuildDEDXTable: Using adapted build
//                               methods of stopping power classes according
//                               to interface change in G4VIonDEDXTable.
//                               Function UpdateCacheValue: Using adapted
//                               ScalingFactorEnergy function according to
//                               interface change in G4VIonDEDXScaling-
//                               Algorithm (AL)
//
// Class description:
//    Ion dE/dx table handler.
//
// Comments:
//
// ===========================================================================

#include <iomanip>

#include "G4IonDEDXHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4VIonDEDXTable.hh"
#include "G4VIonDEDXScalingAlgorithm.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Exp.hh"

//#define PRINT_DEBUG


// #########################################################################

G4IonDEDXHandler::G4IonDEDXHandler(
                            G4VIonDEDXTable* ionTable,
                            G4VIonDEDXScalingAlgorithm* ionAlgorithm,
                            const G4String& name,
                            G4int maxCacheSize,
                            G4bool splines) :
  table(ionTable),
  algorithm(ionAlgorithm),
  tableName(name),
  useSplines(splines),
  maxCacheEntries(maxCacheSize) {

  if(table == 0) {
     G4cerr << "G4IonDEDXHandler::G4IonDEDXHandler() "
            << " Pointer to G4VIonDEDXTable object is null-pointer."
            << G4endl;
  }

  if(algorithm == 0) {
     G4cerr << "G4IonDEDXHandler::G4IonDEDXHandler() "
            << " Pointer to G4VIonDEDXScalingAlgorithm object is null-pointer."
            << G4endl;
  }

  if(maxCacheEntries <= 0) {
     G4cerr << "G4IonDEDXHandler::G4IonDEDXHandler() "
            << " Cache size <=0. Resetting to 5."
            << G4endl;
     maxCacheEntries = 5;
  }
}

// #########################################################################

G4IonDEDXHandler::~G4IonDEDXHandler() {

  ClearCache();

  // All stopping power vectors built according to Bragg's addivitiy rule
  // are deleted. All other stopping power vectors are expected to be
  // deleted by their creator class (sub-class of G4VIonDEDXTable).
  // DEDXTableBraggRule::iterator iter = stoppingPowerTableBragg.begin();
  // DEDXTableBraggRule::iterator iter_end = stoppingPowerTableBragg.end();

  //  for(;iter != iter_end; iter++) delete iter -> second;
  stoppingPowerTableBragg.clear();

  stoppingPowerTable.clear();

  if(table != 0) delete table;
  if(algorithm != 0) delete algorithm;
}

// #########################################################################

G4bool G4IonDEDXHandler::IsApplicable(
                 const G4ParticleDefinition* particle,  // Projectile (ion)
                 const G4Material* material) {          // Target material

  G4bool isApplicable = true;

  if(table == 0 || algorithm == 0) {
     isApplicable = false;
  }
  else {

     G4int atomicNumberIon = particle -> GetAtomicNumber();
     G4int atomicNumberBase =
                algorithm -> AtomicNumberBaseIon(atomicNumberIon, material);

     G4IonKey key = std::make_pair(atomicNumberBase, material);

     DEDXTable::iterator iter = stoppingPowerTable.find(key);
     if(iter == stoppingPowerTable.end()) isApplicable = false;
  }

  return isApplicable;
}

// #########################################################################

G4double G4IonDEDXHandler::GetDEDX(
                 const G4ParticleDefinition* particle,  // Projectile (ion)
                 const G4Material* material,   // Target material
                 G4double kineticEnergy) {     // Kinetic energy of projectile

  G4double dedx = 0.0;

  G4CacheValue value = GetCacheValue(particle, material);

  if(kineticEnergy <= 0.0) dedx = 0.0;
  else if(value.dedxVector != 0) {

     G4bool b;
     G4double factor = value.density;

     factor *= algorithm -> ScalingFactorDEDX(particle,
                                             material,
                                             kineticEnergy);
     G4double scaledKineticEnergy = kineticEnergy * value.energyScaling;

     if(scaledKineticEnergy < value.lowerEnergyEdge) {

        factor *= std::sqrt(scaledKineticEnergy / value.lowerEnergyEdge);
        scaledKineticEnergy = value.lowerEnergyEdge;
     }

     dedx = factor * value.dedxVector -> GetValue(scaledKineticEnergy, b);

     if(dedx < 0.0) dedx = 0.0;
  }
  else dedx = 0.0;

#ifdef PRINT_DEBUG
     G4cout << "G4IonDEDXHandler::GetDEDX() E = "
            << kineticEnergy / MeV << " MeV * "
            << value.energyScaling << " = "
            << kineticEnergy * value.energyScaling / MeV
            << " MeV, dE/dx = " << dedx / MeV * cm << " MeV/cm"
            << ", material = " << material -> GetName()
            << G4endl;
#endif

  return dedx;
}

// #########################################################################

G4bool G4IonDEDXHandler::BuildDEDXTable(
                 const G4ParticleDefinition* particle,  // Projectile (ion)
                 const G4Material* material) {          // Target material

  G4int atomicNumberIon = particle -> GetAtomicNumber();

  G4bool isApplicable = BuildDEDXTable(atomicNumberIon, material);

  return isApplicable;
}


// #########################################################################

G4bool G4IonDEDXHandler::BuildDEDXTable(
                 G4int atomicNumberIon,                // Projectile (ion)
                 const G4Material* material) {         // Target material

  G4bool isApplicable = true;

  if(table == 0 || algorithm == 0) {
     isApplicable = false;
     return isApplicable;
  }

  G4int atomicNumberBase =
                algorithm -> AtomicNumberBaseIon(atomicNumberIon, material);

  // Checking if vector is already built, and returns if this is indeed
  // the case
  G4IonKey key = std::make_pair(atomicNumberBase, material);

  DEDXTable::iterator iter = stoppingPowerTable.find(key);
  if(iter != stoppingPowerTable.end()) return isApplicable;

  // Checking if table contains stopping power vector for given material name
  // or chemical formula
  const G4String& chemFormula = material -> GetChemicalFormula();
  const G4String& materialName = material -> GetName();

  isApplicable = table -> BuildPhysicsVector(atomicNumberBase, chemFormula);

  if(isApplicable) {
     stoppingPowerTable[key] =
              table -> GetPhysicsVector(atomicNumberBase, chemFormula);
     return isApplicable;
  }

  isApplicable = table -> BuildPhysicsVector(atomicNumberBase, materialName);
  if(isApplicable) {
     stoppingPowerTable[key] =
              table -> GetPhysicsVector(atomicNumberBase, materialName);
     return isApplicable;
  }

  // Building the stopping power vector based on Bragg's additivity rule
  const G4ElementVector* elementVector = material -> GetElementVector() ;

  std::vector<G4PhysicsVector*> dEdxTable;

  size_t nmbElements = material -> GetNumberOfElements();

  for(size_t i = 0; i < nmbElements; i++) {

      G4int atomicNumberMat = G4int((*elementVector)[i] -> GetZ());

      isApplicable = table -> BuildPhysicsVector(atomicNumberBase, atomicNumberMat);

      if(isApplicable) {

         G4PhysicsVector* dEdx =
                  table -> GetPhysicsVector(atomicNumberBase, atomicNumberMat);
         dEdxTable.push_back(dEdx);
      }
      else {

         dEdxTable.clear();
         break;
      }
  }

  if(isApplicable) {

     if(dEdxTable.size() > 0) {

        size_t nmbdEdxBins = dEdxTable[0] -> GetVectorLength();
        G4double lowerEdge = dEdxTable[0] -> GetLowEdgeEnergy(0);
        G4double upperEdge = dEdxTable[0] -> GetLowEdgeEnergy(nmbdEdxBins-1);

        G4LPhysicsFreeVector* dEdxBragg =
                    new G4LPhysicsFreeVector(nmbdEdxBins,
                                             lowerEdge,
                                             upperEdge);

        const G4double* massFractionVector = material -> GetFractionVector();

        G4bool b;
        for(size_t j = 0; j < nmbdEdxBins; j++) {

            G4double edge = dEdxTable[0] -> GetLowEdgeEnergy(j);

            G4double value = 0.0;
  	    for(size_t i = 0; i < nmbElements; i++) {

                value += (dEdxTable[i] -> GetValue(edge ,b)) *
                                                       massFractionVector[i];
	    }

            dEdxBragg -> PutValues(j, edge, value);
	}
        dEdxBragg -> SetSpline(useSplines);

#ifdef PRINT_DEBUG
        G4cout << "G4IonDEDXHandler::BuildPhysicsVector() for ion with Z="
               << atomicNumberBase << " in "
               << material -> GetName()
               << G4endl;

        G4cout << *dEdxBragg;
#endif

	stoppingPowerTable[key] = dEdxBragg;
	stoppingPowerTableBragg[key] = dEdxBragg;
     }
  }

  ClearCache();

  return isApplicable;
}

// #########################################################################

G4CacheValue G4IonDEDXHandler::UpdateCacheValue(
              const G4ParticleDefinition* particle,  // Projectile (ion)
              const G4Material* material) {          // Target material

  G4CacheValue value;

  G4int atomicNumberIon = particle -> GetAtomicNumber();
  G4int atomicNumberBase =
                algorithm -> AtomicNumberBaseIon(atomicNumberIon, material);

  G4IonKey key = std::make_pair(atomicNumberBase, material);

  DEDXTable::iterator iter = stoppingPowerTable.find(key);

  if(iter != stoppingPowerTable.end()) {
     value.dedxVector = iter -> second;

     G4double nmbNucleons = G4double(particle -> GetAtomicMass());
     value.energyScaling =
           algorithm -> ScalingFactorEnergy(particle, material) / nmbNucleons;

     size_t nmbdEdxBins = value.dedxVector -> GetVectorLength();
     value.lowerEnergyEdge = value.dedxVector -> GetLowEdgeEnergy(0);

     value.upperEnergyEdge =
                       value.dedxVector -> GetLowEdgeEnergy(nmbdEdxBins-1);
     value.density = material -> GetDensity();
  }
  else {
     value.dedxVector = 0;
     value.energyScaling = 0.0;
     value.lowerEnergyEdge = 0.0;
     value.upperEnergyEdge = 0.0;
     value.density = 0.0;
  }

#ifdef PRINT_DEBUG
  G4cout << "G4IonDEDXHandler::UpdateCacheValue() for "
         << particle -> GetParticleName() << " in "
         << material -> GetName()
         << G4endl;
#endif

  return value;
}

// #########################################################################

G4CacheValue G4IonDEDXHandler::GetCacheValue(
              const G4ParticleDefinition* particle,  // Projectile (ion)
              const G4Material* material) {          // Target material

  G4CacheKey key = std::make_pair(particle, material);

  G4CacheEntry entry;
  CacheEntryList::iterator* pointerIter =
                  (CacheEntryList::iterator*) cacheKeyPointers[key];

  if(!pointerIter) {
      entry.value = UpdateCacheValue(particle, material);

      entry.key = key;
      cacheEntries.push_front(entry);

      CacheEntryList::iterator* pointerIter1 =
                                   new CacheEntryList::iterator();
      *pointerIter1 = cacheEntries.begin();
      cacheKeyPointers[key] = pointerIter1;

      if(G4int(cacheEntries.size()) > maxCacheEntries) {

  	 G4CacheEntry lastEntry = cacheEntries.back();

         void* pointerIter2 = cacheKeyPointers[lastEntry.key];
         CacheEntryList::iterator* listPointerIter =
                      	  (CacheEntryList::iterator*) pointerIter2;

         cacheEntries.erase(*listPointerIter);

         delete listPointerIter;
         cacheKeyPointers.erase(lastEntry.key);
      }
  }
  else {
      entry = *(*pointerIter);
      // Cache entries are currently not re-ordered.
      // Uncomment for activating re-ordering:
      //      cacheEntries.erase(*pointerIter);
      //      cacheEntries.push_front(entry);
      //      *pointerIter = cacheEntries.begin();
  }
  return entry.value;
}

// #########################################################################

void G4IonDEDXHandler::ClearCache() {

  CacheIterPointerMap::iterator iter = cacheKeyPointers.begin();
  CacheIterPointerMap::iterator iter_end = cacheKeyPointers.end();

  for(;iter != iter_end; iter++) {
      void* pointerIter = iter -> second;
      CacheEntryList::iterator* listPointerIter =
                      	  (CacheEntryList::iterator*) pointerIter;

      delete listPointerIter;
  }

  cacheEntries.clear();
  cacheKeyPointers.clear();
}

// #########################################################################

void G4IonDEDXHandler::PrintDEDXTable(
                  const G4ParticleDefinition* particle,  // Projectile (ion)
                  const G4Material* material,  // Target material
                  G4double lowerBoundary,      // Minimum energy per nucleon
                  G4double upperBoundary,      // Maximum energy per nucleon
                  G4int nmbBins,               // Number of bins
                  G4bool logScaleEnergy) {     // Logarithmic scaling of energy

  G4double atomicMassNumber = particle -> GetAtomicMass();
  G4double materialDensity = material -> GetDensity();

  G4cout << "# dE/dx table for " << particle -> GetParticleName()
         << " in material " << material -> GetName()
         << " of density " << materialDensity / g * cm3
         << " g/cm3"
         << G4endl
         << "# Projectile mass number A1 = " << atomicMassNumber
         << G4endl
         << "# Energy range (per nucleon) of tabulation: "
         << GetLowerEnergyEdge(particle, material) / atomicMassNumber / MeV
         << " - "
         << GetUpperEnergyEdge(particle, material) / atomicMassNumber / MeV
         << " MeV"
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;
  G4cout << "#"
         << std::setw(13) << std::right << "E"
         << std::setw(14) << "E/A1"
         << std::setw(14) << "dE/dx"
         << std::setw(14) << "1/rho*dE/dx"
         << G4endl;
  G4cout << "#"
         << std::setw(13) << std::right << "(MeV)"
         << std::setw(14) << "(MeV)"
         << std::setw(14) << "(MeV/cm)"
         << std::setw(14) << "(MeV*cm2/mg)"
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;

  //G4CacheValue value = GetCacheValue(particle, material);

  G4double energyLowerBoundary = lowerBoundary * atomicMassNumber;
  G4double energyUpperBoundary = upperBoundary * atomicMassNumber;

  if(logScaleEnergy) {

     energyLowerBoundary = std::log(energyLowerBoundary);
     energyUpperBoundary = std::log(energyUpperBoundary);
  }

  G4double deltaEnergy = (energyUpperBoundary - energyLowerBoundary) /
                                                           G4double(nmbBins);

  G4cout.precision(6);
  for(int i = 0; i < nmbBins + 1; i++) {

      G4double energy = energyLowerBoundary + i * deltaEnergy;
      if(logScaleEnergy) energy = G4Exp(energy);

      G4double loss = GetDEDX(particle, material, energy);

      G4cout << std::setw(14) << std::right << energy / MeV
             << std::setw(14) << energy / atomicMassNumber / MeV
             << std::setw(14) << loss / MeV * cm
             << std::setw(14) << loss / materialDensity / (MeV*cm2/(0.001*g))
             << G4endl;
  }
}

// #########################################################################

G4double G4IonDEDXHandler::GetLowerEnergyEdge(
                 const G4ParticleDefinition* particle,  // Projectile (ion)
                 const G4Material* material) {          // Target material

  G4double edge = 0.0;

  G4CacheValue value = GetCacheValue(particle, material);

  if(value.energyScaling > 0)
          edge = value.lowerEnergyEdge / value.energyScaling;

  return edge;
}

// #########################################################################

G4double G4IonDEDXHandler::GetUpperEnergyEdge(
                 const G4ParticleDefinition* particle,  // Projectile (ion)
                 const G4Material* material) {          // Target material

  G4double edge = 0.0;

  G4CacheValue value = GetCacheValue(particle, material);

  if(value.energyScaling > 0)
          edge = value.upperEnergyEdge / value.energyScaling;

  return edge;
}

// #########################################################################

G4String G4IonDEDXHandler::GetName() {

  return tableName;
}

// #########################################################################
