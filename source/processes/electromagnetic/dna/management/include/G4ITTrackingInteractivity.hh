#ifndef G4ITTRACKINGINTERACTIVITY_HH
#define G4ITTRACKINGINTERACTIVITY_HH

#include "G4String.hh"

class G4Track;
class G4Step;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4ITTrackingInteractivity
{
protected :
    int fVerboseLevel;

public:
    G4ITTrackingInteractivity(){;}
    virtual ~G4ITTrackingInteractivity(){;}

    virtual void StartTracking(G4Track*){;}
    virtual void AppendStep(G4Track* /*track*/, G4Step* /*step*/){;}
    virtual void EndTracking(G4Track*){;}

    virtual void TrackBanner(G4Track* /*track*/, const G4String& message = "");
    inline void SetVerbose(int);
};

inline void G4ITTrackingInteractivity::SetVerbose(int flag)
{
    fVerboseLevel = flag;
}

#endif // G4ITTRACKINGINTERACTIVITY_HH
