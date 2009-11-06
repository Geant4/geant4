/////////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointStackingAction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	 -April  2008 First implementation  by L. Desorgher 
// 		 -4-11-2009 Adding the possibility to use user adjoint stacking  action,  L. Desorgher
//	 	 
//
//-------------------------------------------------------------
//	Documentation:
//		Stacking action used in the adjoint simulation. It is responsible to kill a primary forward particle before it is traked  in the forwrad phase 
//		if the last adjoint particle did not reach the adjoint surface. Was needed for the new design where the G4AdjointSimManager is no more an extension
//		of the G4RunManager. If the primary particles is not killed before being tracked in the sensitive geometry, the User Stacking action 
//		can be used duiring the forward phase if specified by the method  G4AdjointSimManager::UseUserStackingAction(Bool). 
//	
//
//
//	
//
//		
//
#ifndef G4AdjointStackingAction_h
#define G4AdjointStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"



class G4Track;
class G4ParticleDefinition;



class G4AdjointStackingAction : public G4UserStackingAction
{
  public:
    G4AdjointStackingAction();
    virtual ~G4AdjointStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();
    inline void SetUserFwdStackingAction(G4UserStackingAction* anAction){theFwdStackingAction = anAction;}
    inline void SetUserAdjointStackingAction(G4UserStackingAction* anAction){theUserAdjointStackingAction = anAction;}
    
    inline void SetKillTracks(G4bool aBool){kill_tracks =aBool;}
    inline void SetAdjointMode(G4bool aBool){adjoint_mode=aBool;}
  
  private:
    G4UserStackingAction* theFwdStackingAction;
    G4UserStackingAction* theUserAdjointStackingAction;
    G4bool kill_tracks, adjoint_mode;
  
    
};
  

#endif

