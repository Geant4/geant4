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

#include "G4INCLStandaloneParticleDataSource.hh"

namespace G4INCL {
  StandaloneParticleDataSource::StandaloneParticleDataSource()
  {}

  StandaloneParticleDataSource::~StandaloneParticleDataSource()
  {}

  std::string StandaloneParticleDataSource::getPDSName(ParticleType /*t*/)
  {
    return std::string("NAMES_NOT_IMPLEMENTED_YET");
  }

  std::string StandaloneParticleDataSource::getPDSName(G4int /*A*/, G4int /*Z*/)
  {
    return std::string("NAMES_NOT_IMPLEMENTED_YET");
  }

  G4double StandaloneParticleDataSource::getMass(ParticleType t)
  {
    const G4double mN = 938.2796;
    const G4double mPi = 138.0;
    if(t == Proton || t == Neutron) {
      return mN;
    } else if(t == PiPlus || t == PiMinus || t == PiZero) {
      return mPi;
    } else {
      ERROR("Standalone PDS: No mass for unknown particle." << std::endl);
      return 0.0;
    }
  }

  G4double StandaloneParticleDataSource::getMass(G4int A, G4int Z)
  {
    return (A - Z) * getMass(Neutron) + Z * getMass(Proton);
  }
}
