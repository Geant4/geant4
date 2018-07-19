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
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLCDPP_hh
#define G4INCLCDPP_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPauli.hh"
#include "G4INCLINuclearPotential.hh"

namespace G4INCL {
  class CDPP : public IPauli {
  public:
    CDPP();
    ~CDPP();

    G4bool isBlocked(ParticleList const &, Nucleus const * const);

    inline void processOneParticle(Particle const * const p) {
      if(p->isBaryon()) {
        const G4double Tf = thePotential->getFermiEnergy(p);
        const G4double T = p->getKineticEnergy();

        if(T > Tf) {
          const G4double sep = thePotential->getSeparationEnergy(p);
          Sk += sep;
        } else {
          TbelowTf += T - p->getPotentialEnergy();
        }
      } else if(p->isMeson() || p->isResonance()) {
        const G4double sep = thePotential->getSeparationEnergy(p);
        Sk += sep;
      }
    }

  protected:
    G4double Sk;
    G4double TbelowTf;
    NuclearPotential::INuclearPotential const * thePotential;

  };
}

#endif
