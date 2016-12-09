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
// $Id: G4AdjointStackingAction.hh 98744 2016-08-09 13:21:26Z gcosmo $
//
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
class G4AdjointTrackingAction;

class G4AdjointStackingAction : public G4UserStackingAction
{
  public:
    G4AdjointStackingAction(G4AdjointTrackingAction* anAction);
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
    G4bool reclassification_stage,first_reclassification_stage,kill_tracks,adjoint_mode;
    G4AdjointTrackingAction* theAdjointTrackingAction;
};

#endif

