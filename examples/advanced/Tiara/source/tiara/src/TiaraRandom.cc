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


