// $Id: G4MCTGenParticle.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// ====================================================================
//
//   G4MCTGenParticle.hh
//
// ====================================================================
#ifndef MCT_GEN_PARTICLE_H
#define MCT_GEN_PARTICLE_H

#include <utility>
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenParticle.h"

class HepMC::GenEvent;
class HepMC::GenParticle;

typedef std::pair<HepMC::GenEvent*, HepMC::GenParticle*> G4MCTGenParticle;

#endif

