//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
 * G4VTrackState.hh
 *
 *  Created on: 8 avr. 2013
 *      Author: kara
 */

#ifndef G4VTRACKSTATE_HH_
#define G4VTRACKSTATE_HH_

#include <CLHEP/Utility/memory.h>

class G4VTrackStateID
{
public :
	virtual int GetStateID() = 0;
protected :
	G4VTrackStateID(){;}
	virtual ~G4VTrackStateID(){;}

	static int fgLastID;
};

class G4TrackState;

template <class T>
class G4TrackStateID : public G4VTrackStateID
{
public :
	static CLHEP::shared_ptr<G4VTrackStateID> Get()
	{
		if(!fgTrackStateID)
		{
			new G4TrackStateID<T>;
		}

		return fgTrackStateID;
	}

	~G4TrackStateID<T>(){;}

	virtual int GetStateID()
	{
		return fID;
	}

	static int GetID()
	{
		if(fgTrackStateID)
			return fgTrackStateID->fID;
		else
			return -1; // TODO Exception
	}

private:
	G4TrackStateID<T>() : fID(fgLastID)
	{
		fgLastID++;
		fgTrackStateID = this;
	}

	static CLHEP::shared_ptr<G4TrackStateID<T> > fgTrackStateID;
	const int fID;
};


class G4TrackState
{
public:
	template <class T> static G4TrackState* Create()
	{
		G4TrackState* output = new G4TrackState(G4TrackStateID<T>::Get());
		return output;
	}

	int GetID()
	{
		return fpTrackStateID->GetStateID();
	}

protected:
	virtual ~G4TrackState(){}

private :
	G4TrackState(CLHEP::shared_ptr<G4VTrackStateID> _trackStateID) : fpTrackStateID(_trackStateID)
	{}

	CLHEP::shared_ptr<G4VTrackStateID> fpTrackStateID;
};

#endif /* G4VTRACKSTATE_HH_ */
