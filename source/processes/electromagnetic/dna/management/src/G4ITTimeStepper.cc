/*
 * G4ITTimeStepper.cpp
 *
 *  Created on: 4 juin 2014
 *      Author: kara
 */

#include "G4ITTimeStepper.hh"
#include "G4ITStepManager.hh"

G4ThreadLocal G4ITTimeStepper* G4ITTimeStepper::fpInstance;

G4ITTimeStepper::G4ITTimeStepper() {
	fpInstance = this;
}

G4ITTimeStepper::~G4ITTimeStepper() {
	fpInstance = 0;
}

G4ITTimeStepper* G4ITTimeStepper::Instance(){
 return G4ITStepManager::Instance();
}
