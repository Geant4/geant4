#ifndef CML2TrackingActionH
#define CML2TrackingActionH

#include "globals.hh"
#include "G4UserTrackingAction.hh"

class G4Navigator;
class G4Track; 

class CML2TrackingAction  :  public G4UserTrackingAction 
{
public:
	CML2TrackingAction(void);
	virtual ~CML2TrackingAction(void);
     virtual void PreUserTrackingAction(const G4Track*);

private:
      G4Navigator* gNavigator; 
};

#endif
