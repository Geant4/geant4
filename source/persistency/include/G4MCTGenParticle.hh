// $Id: G4MCTGenParticle.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// ====================================================================
//
//   G4MCTGenParticle.hh
//
// ====================================================================
#ifndef MCT_GEN_PARTICLE_H
#define MCT_GEN_PARTICLE_H

#include "G4Types.hh"
#include "g4std/algorithm"
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenParticle.h"

class HepMC::GenEvent;
class HepMC::GenParticle;

typedef G4std::pair<HepMC::GenEvent*, HepMC::GenParticle*> G4MCTGenParticle;

#endif

