/*
 * G4MolecularGun.cc
 *
 *  Created on: 28 janv. 2014
 *      Author: kara
 */

#include "G4ITGun.hh"
#include "G4ITTrackHolder.hh"

G4ITGun::G4ITGun() {
	// TODO Auto-generated constructor stub

}

G4ITGun::~G4ITGun() {
	// TODO Auto-generated destructor stub
}

void G4ITGun::PushTrack(G4Track* track)
{
	G4VITTrackHolder::Instance()->Push(track);
}
