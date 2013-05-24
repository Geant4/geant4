/*
 * G4VTrackState.cc
 *
 *  Created on: 8 avr. 2013
 *      Author: kara
 */

#include "G4VTrackState.hh"

int G4VTrackStateID::fgLastID = 0;

template <class T> CLHEP::shared_ptr<G4TrackStateID<T> > G4TrackStateID<T>::fgTrackStateID (0);
