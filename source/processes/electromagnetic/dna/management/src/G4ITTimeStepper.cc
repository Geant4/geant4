/*
 * G4ITTimeStepper.cpp
 *
 *  Created on: 4 juin 2014
 *      Author: kara
 */

#include "G4ITTimeStepper.hh"

G4ThreadLocal G4ITTimeStepper* G4ITTimeStepper::fpInstance;

G4ITTimeStepper::G4ITTimeStepper() {
	fpInstance = this;
}

G4ITTimeStepper::~G4ITTimeStepper() {
	fpInstance = 0;
}

