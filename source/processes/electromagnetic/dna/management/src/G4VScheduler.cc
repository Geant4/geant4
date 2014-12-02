/*
 * G4ITTimeStepper.cpp
 *
 *  Created on: 4 juin 2014
 *      Author: kara
 */

#include <G4Scheduler.hh>
#include "G4VScheduler.hh"

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
  //========================================================================
  // For now G4VScheduler::Instance() returns G4ITScheduler
  //========================================================================
  return G4Scheduler::Instance();
}

void G4VScheduler::Process()
{
  G4Exception("G4VScheduler::Process()",
              "CONCRETE_OBJECT_MISSING",
              FatalException,
              "A concrete implementation of G4VScheduler is expected to be "
              "instantiate before starting the run, have you correctly "
              "initialized the chemistry ?");
}
