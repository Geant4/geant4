///////////////////////
//G4RTTrackingAction.hh
///////////////////////


#ifndef G4RTTrackingAction_h
#define G4RTTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class G4Track;

///////////////////////////
class G4RTTrackingAction : public G4UserTrackingAction
///////////////////////////
{

//--------
   public:
//--------

// Constructor & Destructor
   G4RTTrackingAction(){;}
   virtual ~G4RTTrackingAction(){;}

// Member functions
   virtual void PreUserTrackingAction(const G4Track* aTrack);
   virtual void PostUserTrackingAction(const G4Track* aTrack){;}


};

#endif


