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
// INCL++ revision: v5.1.1
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNuclearDensityFactory.hh"

namespace G4INCL {

  std::map<G4int,NuclearDensity*> NuclearDensityFactory::nuclearDensityCache;

  NuclearDensity* NuclearDensityFactory::createDensity(const G4int A, const G4int Z, const G4bool hardFermiSphere/*=true*/) {
    const G4int nuclideID = (1000*Z + A)*(hardFermiSphere?(-1):1); // MCNP-style nuclide IDs
    const std::map<G4int,NuclearDensity*>::const_iterator mapEntry = nuclearDensityCache.find(nuclideID);
    if(mapEntry == nuclearDensityCache.end()) {
      IFunction1D *densityFunction = NuclearDensityFactory::createDensityFunction(A, Z);
      if(!densityFunction)
        return NULL;
      NuclearDensity *density = new NuclearDensity(A, Z, densityFunction, hardFermiSphere);
      nuclearDensityCache[nuclideID] = density;
      return density;
    } else {
      return mapEntry->second;
    }
  }

  IFunction1D* NuclearDensityFactory::createDensityFunction(const G4int A, const G4int Z) {
    G4double radius = ParticleTable::getNuclearRadius(A, Z);
    G4double diffuseness = ParticleTable::getSurfaceDiffuseness(A, Z);
    G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(A, Z);

    if(A > 19) {
      return new DerivWoodsSaxon(radius, maximumRadius, diffuseness);
    } else if(A <= 19 && A > 6) {
      return new DerivModifiedHarmonicOscillator(radius, maximumRadius, diffuseness);
    } else if(A >= 2 && A <= 6) { // Gaussian distribution for light nuclei
      return new DerivGaussian(radius, maximumRadius, diffuseness);
    } else {
      ERROR("No nuclear density function for target A = "
        << A << " Z = " << Z << std::endl);
    }
    return 0;
  }

}
