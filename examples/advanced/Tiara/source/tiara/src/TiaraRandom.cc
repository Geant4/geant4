// $Id: TiaraRandom.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraRandom.hh"

#include "Randomize.hh"

void setRandomSeed(long seed) {
  HepRandom::setTheSeed(seed);
}

void setRandomStatus(const char *filename){
    HepRandom::restoreEngineStatus(filename);
}

void saveRandomStatus(const char *filename){
  HepRandom::saveEngineStatus(filename);
}


