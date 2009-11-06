#include "G4AdjointStackingAction.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4ios.hh"


G4AdjointStackingAction::G4AdjointStackingAction()
{ theFwdStackingAction =0;
  theUserAdjointStackingAction =0;
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointStackingAction::~G4AdjointStackingAction()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4ClassificationOfNewTrack  G4AdjointStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{  
   G4ClassificationOfNewTrack classification = fUrgent;
   if ( kill_tracks) classification=fKill;
   else if (!adjoint_mode && theFwdStackingAction)   classification =  theFwdStackingAction->ClassifyNewTrack(aTrack); 
   else if (adjoint_mode && theUserAdjointStackingAction)   classification =  theUserAdjointStackingAction->ClassifyNewTrack(aTrack);
   return classification;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointStackingAction::NewStage()
{ if ( !adjoint_mode && theFwdStackingAction)  theFwdStackingAction->NewStage(); 
  else if (adjoint_mode && theUserAdjointStackingAction)   theUserAdjointStackingAction->NewStage();
}
////////////////////////////////////////////////////////////////////////////////
//    
void G4AdjointStackingAction::PrepareNewEvent()
{ if ( !adjoint_mode && theFwdStackingAction)  theFwdStackingAction->PrepareNewEvent();
  else if (adjoint_mode && theUserAdjointStackingAction)   theUserAdjointStackingAction->PrepareNewEvent();
}



