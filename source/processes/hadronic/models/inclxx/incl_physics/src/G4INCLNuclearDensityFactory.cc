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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLNDFWoodsSaxon.hh"
#include "G4INCLNDFModifiedHarmonicOscillator.hh"
#include "G4INCLNDFGaussian.hh"
#include "G4INCLNDFParis.hh"
#include "G4INCLNDFHardSphere.hh"

namespace G4INCL {

  namespace NuclearDensityFactory {

    namespace {

      G4ThreadLocal std::map<G4int,NuclearDensity const *> *nuclearDensityCache = NULL;
      G4ThreadLocal std::map<G4int,InverseInterpolationTable*> *rpCorrelationTableCache = NULL;
      G4ThreadLocal std::map<G4int,InverseInterpolationTable*> *rCDFTableCache = NULL;
      G4ThreadLocal std::map<G4int,InverseInterpolationTable*> *pCDFTableCache = NULL;

    }

    NuclearDensity const *createDensity(const G4int A, const G4int Z) {
      if(!nuclearDensityCache)
        nuclearDensityCache = new std::map<G4int,NuclearDensity const *>;

      const G4int nuclideID = 1000*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,NuclearDensity const *>::const_iterator mapEntry = nuclearDensityCache->find(nuclideID);
      if(mapEntry == nuclearDensityCache->end()) {
        InverseInterpolationTable *rpCorrelationTableProton = createRPCorrelationTable(Proton, A, Z);
        InverseInterpolationTable *rpCorrelationTableNeutron = createRPCorrelationTable(Neutron, A, Z);
        if(!rpCorrelationTableProton || !rpCorrelationTableNeutron)
          return NULL;
        NuclearDensity const *density = new NuclearDensity(A, Z, rpCorrelationTableProton, rpCorrelationTableNeutron);
        (*nuclearDensityCache)[nuclideID] = density;
        return density;
      } else {
        return mapEntry->second;
      }
    }

    InverseInterpolationTable *createRPCorrelationTable(const ParticleType t, const G4int A, const G4int Z) {
// assert(t==Proton || t==Neutron);

      if(!rpCorrelationTableCache)
        rpCorrelationTableCache = new std::map<G4int,InverseInterpolationTable*>;

      const G4int nuclideID = ((t==Proton) ? 1000 : -1000)*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,InverseInterpolationTable*>::const_iterator mapEntry = rpCorrelationTableCache->find(nuclideID);
      if(mapEntry == rpCorrelationTableCache->end()) {

        IFunction1D *rpCorrelationFunction;
        if(A > 19) {
          G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(t, A, Z);
          const G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rpCorrelationFunction = new NuclearDensityFunctions::WoodsSaxonRP(radius, maximumRadius, diffuseness);
        } else if(A <= 19 && A > 6) {
          const G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          const G4double diffuseness = ParticleTable::getSurfaceDiffuseness(t, A, Z);
          const G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rpCorrelationFunction = new NuclearDensityFunctions::ModifiedHarmonicOscillatorRP(radius, maximumRadius, diffuseness);
        } else if(A <= 6 && A > 1) { // Gaussian distribution for light nuclei
          const G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          const G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rpCorrelationFunction = new NuclearDensityFunctions::GaussianRP(maximumRadius, Math::oneOverSqrtThree * radius);
        } else {
          INCL_ERROR("No r-p correlation function for " << ((t==Proton) ? "protons" : "neutrons") << " in A = "
                << A << " Z = " << Z << std::endl);
          return NULL;

        }

        class InverseCDFOneThird : public IFunction1D {
          public:
            InverseCDFOneThird(IFunction1D const * const f) :
              IFunction1D(f->getXMinimum(), f->getXMaximum()),
              theFunction(f),
              normalisation(1./theFunction->integrate(xMin,xMax))
          {}

            G4double operator()(const G4double x) const {
              return Math::pow13(normalisation * theFunction->integrate(xMin,x));
            }
          private:
            IFunction1D const * const theFunction;
            const G4double normalisation;
        } *theInverseCDFOneThird = new InverseCDFOneThird(rpCorrelationFunction);

        InverseInterpolationTable *theTable = new InverseInterpolationTable(*theInverseCDFOneThird);
        delete theInverseCDFOneThird;
        delete rpCorrelationFunction;
        INCL_DEBUG("Creating r-p correlation function for " << ((t==Proton) ? "protons" : "neutrons") << " in A=" << A << ", Z=" << Z << ":"
              << std::endl << theTable->print() << std::endl);

        (*rpCorrelationTableCache)[nuclideID] = theTable;
        return theTable;
      } else {
        return mapEntry->second;
      }
    }

    InverseInterpolationTable *createRCDFTable(const ParticleType t, const G4int A, const G4int Z) {
// assert(t==Proton || t==Neutron);

      if(!rCDFTableCache)
        rCDFTableCache = new std::map<G4int,InverseInterpolationTable*>;

      const G4int nuclideID = ((t==Proton) ? 1000 : -1000)*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,InverseInterpolationTable*>::const_iterator mapEntry = rCDFTableCache->find(nuclideID);
      if(mapEntry == rCDFTableCache->end()) {

        IFunction1D *rDensityFunction;
        if(A > 19) {
          G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(t, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rDensityFunction = new NuclearDensityFunctions::WoodsSaxon(radius, maximumRadius, diffuseness);
        } else if(A <= 19 && A > 6) {
          G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(t, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rDensityFunction = new NuclearDensityFunctions::ModifiedHarmonicOscillator(radius, maximumRadius, diffuseness);
        } else if(A <= 6 && A > 2) { // Gaussian distribution for light nuclei
          G4double radius = ParticleTable::getRadiusParameter(t, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(t, A, Z);
          rDensityFunction = new NuclearDensityFunctions::Gaussian(maximumRadius, Math::oneOverSqrtThree * radius);
        } else if(A == 2 && Z == 1) { // density from the Paris potential for deuterons
          rDensityFunction = new NuclearDensityFunctions::ParisR();
        } else {
          INCL_ERROR("No nuclear density function for target A = "
                << A << " Z = " << Z << std::endl);
          return NULL;
        }

        InverseInterpolationTable *theTable = rDensityFunction->inverseCDFTable();
        delete rDensityFunction;
        INCL_DEBUG("Creating inverse position CDF for A=" << A << ", Z=" << Z << ":" <<
              std::endl << theTable->print() << std::endl);

        (*rCDFTableCache)[nuclideID] = theTable;
        return theTable;
      } else {
        return mapEntry->second;
      }
    }

    InverseInterpolationTable *createPCDFTable(const ParticleType t, const G4int A, const G4int Z) {
// assert(t==Proton || t==Neutron);

      if(!pCDFTableCache)
        pCDFTableCache = new std::map<G4int,InverseInterpolationTable*>;

      const G4int nuclideID = ((t==Proton) ? 1000 : -1000)*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,InverseInterpolationTable*>::const_iterator mapEntry = pCDFTableCache->find(nuclideID);
      if(mapEntry == pCDFTableCache->end()) {
        IFunction1D *pDensityFunction;
        if(A > 19) {
          const G4double theFermiMomentum = ParticleTable::getFermiMomentum(A, Z);
          pDensityFunction = new NuclearDensityFunctions::HardSphere(theFermiMomentum);
        } else if(A <= 19 && A > 2) { // Gaussian distribution for light nuclei
          const G4double momentumRMS = Math::oneOverSqrtThree * ParticleTable::getMomentumRMS(A, Z);
          pDensityFunction = new NuclearDensityFunctions::Gaussian(5.*momentumRMS, momentumRMS);
        } else if(A == 2 && Z == 1) { // density from the Paris potential for deuterons
          pDensityFunction = new NuclearDensityFunctions::ParisP();
        } else {
          INCL_ERROR("No nuclear density function for target A = "
                << A << " Z = " << Z << std::endl);
          return NULL;
        }

        InverseInterpolationTable *theTable = pDensityFunction->inverseCDFTable();
        delete pDensityFunction;
        INCL_DEBUG("Creating inverse momentum CDF for A=" << A << ", Z=" << Z << ":" <<
              std::endl << theTable->print() << std::endl);

        (*pCDFTableCache)[nuclideID] = theTable;
        return theTable;
      } else {
        return mapEntry->second;
      }
    }

    void addRPCorrelationToCache(const G4int A, const G4int Z, const ParticleType t, InverseInterpolationTable * const table) {
// assert(t==Proton || t==Neutron);

      if(!rpCorrelationTableCache)
        rpCorrelationTableCache = new std::map<G4int,InverseInterpolationTable*>;

      const G4int nuclideID = ((t==Proton) ? 1000 : -1000)*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,InverseInterpolationTable*>::const_iterator mapEntry = rpCorrelationTableCache->find(nuclideID);
      if(mapEntry != rpCorrelationTableCache->end())
        delete mapEntry->second;

      (*rpCorrelationTableCache)[nuclideID] = table;
    }

    void addDensityToCache(const G4int A, const G4int Z, NuclearDensity * const density) {
      if(!nuclearDensityCache)
        nuclearDensityCache = new std::map<G4int,NuclearDensity const *>;

      const G4int nuclideID = 1000*Z + A; // MCNP-style nuclide IDs
      const std::map<G4int,NuclearDensity const *>::const_iterator mapEntry = nuclearDensityCache->find(nuclideID);
      if(mapEntry != nuclearDensityCache->end())
        delete mapEntry->second;

      (*nuclearDensityCache)[nuclideID] = density;
    }

    void clearCache() {

      if(nuclearDensityCache) {
        for(std::map<G4int,NuclearDensity const *>::const_iterator i = nuclearDensityCache->begin(); i!=nuclearDensityCache->end(); ++i)
          delete i->second;
        nuclearDensityCache->clear();
        delete nuclearDensityCache;
        nuclearDensityCache = NULL;
      }

      if(rpCorrelationTableCache) {
        for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = rpCorrelationTableCache->begin(); i!=rpCorrelationTableCache->end(); ++i)
          delete i->second;
        rpCorrelationTableCache->clear();
        delete rpCorrelationTableCache;
        rpCorrelationTableCache = NULL;
      }

      if(rCDFTableCache) {
        for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = rCDFTableCache->begin(); i!=rCDFTableCache->end(); ++i)
          delete i->second;
        rCDFTableCache->clear();
        delete rCDFTableCache;
        rCDFTableCache = NULL;
      }

      if(pCDFTableCache) {
        for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = pCDFTableCache->begin(); i!=pCDFTableCache->end(); ++i)
          delete i->second;
        pCDFTableCache->clear();
        delete pCDFTableCache;
        pCDFTableCache = NULL;
      }
    }

  } // namespace NuclearDensityFactory

} // namespace G4INCL
