/*
 * G4VTrackState.hh
 *
 *  Created on: 8 avr. 2013
 *      Author: kara
 */

#ifndef G4VTRACKSTATE_HH_
#define G4VTRACKSTATE_HH_

class G4TrackStateCounter
{
protected :
	static int fgLastID;
};

template <class T>
class G4TrackStateID : public G4TrackStateCounter
{
public :
	static int Create()
	{
		if(fgTrackStateID)
		{
			//			G4Exception();
		}
		else
		{
			new G4TrackStateID<T>;
		}

		return fgTrackStateID->fID;
	}
	~G4TrackStateID<T>(){;}
	static int GetID() { return fID; }

protected:
	G4TrackStateID<T>() : fID(fgLastID)
	{
		fgLastID++;
		fgTrackStateID = this;
	}
	G4TrackStateID<T>* fgTrackStateID;
	static const int fID;
};


class G4VTrackState
{
public:
	virtual ~G4VTrackState();

protected:
	G4VTrackState();
	virtual int GetID() = 0;
};

#endif /* G4VTRACKSTATE_HH_ */
