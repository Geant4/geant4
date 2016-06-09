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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNuclearDensityFactory.hh"

namespace G4INCL {

  NuclearDensity* NuclearDensityFactory::createDensity(G4int A, G4int Z) {
    IFunction1D *densityFunction = NuclearDensityFactory::createDensityFunction(A, Z);
    NuclearDensity *density = new NuclearDensity(A, Z, densityFunction);
    return density;
  }

  IFunction1D* NuclearDensityFactory::createDensityFunction(G4int A, G4int Z) {
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

  // We will not construct any instances of this class
  NuclearDensityFactory::NuclearDensityFactory() {}
  NuclearDensityFactory::~NuclearDensityFactory() {}

}
