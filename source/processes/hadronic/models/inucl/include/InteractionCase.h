#ifndef INTERACTION_CASE_H
#define INTERACTION_CASE_H

#include "InuclParticle.h"
#include "pair.h"

class InteractionCase {

public:

InteractionCase() { 
  bultag = pair<InuclParticle*,InuclParticle*>(0,0);
  inter_case = 0;
};

InteractionCase(InuclParticle* part1, InuclParticle* part2, int ic) {
  setBulletTarget(part1,part2);
  setInterCase(ic);
}; 

void setBulletTarget(InuclParticle* part1, InuclParticle* part2) {
  bultag = pair<InuclParticle*,InuclParticle*>(part1,part2);
};

void setInterCase(int ic) { inter_case = ic; };

InuclParticle* getBullet() const { return bultag.first; };

InuclParticle* getTarget() const { return bultag.second; };

int getInterCase() const { return inter_case; };

private:

pair<InuclParticle*,InuclParticle*> bultag;

int inter_case;

};

#endif // INTERACTION_CASE_H 
