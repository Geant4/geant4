 #ifndef Em10TrackingAction_h
 #define Em10TrackingAction_h

 #include "G4UserTrackingAction.hh"

 #include "globals.hh"

 class Em10TrackingAction : public G4UserTrackingAction {

  public:
    Em10TrackingAction();
    virtual ~Em10TrackingAction(){};

    virtual void PreUserTrackingAction(const G4Track*);

 };

 #endif

