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

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLICoulomb.hh"
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


#ifndef G4INCLPbarAtrestEntry_hh
#define G4INCLPbarAtrestEntry_hh 1

namespace G4INCL {

  class FinalState;

  class PbarAtrestEntryChannel : public IChannel {
  public:
    PbarAtrestEntryChannel(Nucleus *n, Particle *p);
    virtual ~PbarAtrestEntryChannel();

    void fillFinalState(FinalState *fs);

    ParticleList makeMesonStar();
    G4double PbarCoulombicCascadeEnergy(G4int A, G4int Z);
    G4double n_annihilation(G4int A, G4int Z);
    IAvatarList bringMesonStar(ParticleList const &pL, Nucleus * const n);
    G4bool ProtonIsTheVictim();
    ThreeVector getAnnihilationPosition();
  
    G4double fctrl(const G4double arg1);
    G4double r1(const G4int n);
    G4double r2(const G4int n);
    G4double r3(G4double x, const G4int n);
    G4double r4(G4double x, const G4int n);
    G4double radial_wavefunction(G4double x, const G4int n);
    G4double densityP(G4double x);
    G4double densityN(G4double x);
    G4double overlapP(G4double &x, const G4int n);
    G4double overlapN(G4double &x, const G4int n);
    G4double read_file(std::string filename, std::vector<G4double>& probabilities, std::vector<std::vector<std::string>>& particle_types);
    G4int findStringNumber(G4double rdm, std::vector<G4double> yields);
    
  private:
    Nucleus *theNucleus;
    Particle *theParticle;


    INCL_DECLARE_ALLOCATION_POOL(PbarAtrestEntryChannel)
  };
}

#endif
