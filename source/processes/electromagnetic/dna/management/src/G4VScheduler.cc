/*
 * G4ITTimeStepper.cpp
 *
 *  Created on: 4 juin 2014
 *      Author: kara
 */

#include "G4VScheduler.hh"
#include "G4ITScheduler.hh"

G4ThreadLocal G4VScheduler* G4VScheduler::fpInstance;

G4VScheduler::G4VScheduler()
{
  fpInstance = this;
}

G4VScheduler::~G4VScheduler()
{
  fpInstance = 0;
}

G4VScheduler* G4VScheduler::Instance()
{
  return G4ITScheduler::Instance();
}
