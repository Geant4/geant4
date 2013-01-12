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

#ifndef G4INCLNuclearDensityFactory_hh
#define G4INCLNuclearDensityFactory_hh 1

#include "G4INCLNuclearDensity.hh"
#include "G4INCLParticleSampler.hh"
#include "G4INCLParticleTable.hh"
#include <map>

namespace G4INCL {

  class NuclearDensityFactory {
  public:
    static NuclearDensity *createDensity(const G4int A, const G4int Z);

    static InverseInterpolationTable *createRPCorrelationTable(const G4int A, const G4int Z);

    static InverseInterpolationTable *createRCDFTable(const G4int A, const G4int Z);

    static InverseInterpolationTable *createPCDFTable(const G4int A, const G4int Z);

    static ParticleSampler *createParticleSampler(const G4int A, const G4int Z);

    static void clearCache() {
      for(std::map<G4int,NuclearDensity*>::const_iterator i = nuclearDensityCache.begin(); i!=nuclearDensityCache.end(); ++i)
        delete i->second;
      nuclearDensityCache.clear();

      for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = rpCorrelationTableCache.begin(); i!=rpCorrelationTableCache.end(); ++i)
        delete i->second;
      rpCorrelationTableCache.clear();

      for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = rCDFTableCache.begin(); i!=rCDFTableCache.end(); ++i)
        delete i->second;
      rCDFTableCache.clear();

      for(std::map<G4int,InverseInterpolationTable*>::const_iterator i = pCDFTableCache.begin(); i!=pCDFTableCache.end(); ++i)
        delete i->second;
      pCDFTableCache.clear();
    }

  protected:
    // We will not construct any instances of this class
    NuclearDensityFactory() {}
    ~NuclearDensityFactory() {}

    static std::map<G4int,NuclearDensity*> nuclearDensityCache;

    static std::map<G4int,InverseInterpolationTable*> rpCorrelationTableCache;
    static std::map<G4int,InverseInterpolationTable*> rCDFTableCache;
    static std::map<G4int,InverseInterpolationTable*> pCDFTableCache;

  };
}

#endif
