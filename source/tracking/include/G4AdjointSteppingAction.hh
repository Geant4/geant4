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
// $Id: G4AdjointSteppingAction.hh 102436 2017-01-27 08:28:53Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointSteppingAction
//	Author:       	L. Desorgher
//	Contract:	ESA contract 21435/08/NL/AT
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//		 -15/01/2007 Creation by  L. Desorgher 
// 		 -01/11/2009 Some cleaning and adding of documentation  for the first Release in the Geant4 toolkit, L. Desorgher 
//		 -04/11/2009 Adding the possibility to use user stepping action,  L. Desorgher
//		 -20/11/2009 Correct the stop of adjoint particles tracking when it reenters the adjoint source.
//		 
//
//-------------------------------------------------------------
//	Documentation:
//		Stepping action used in the adjoint simulation. 
//		It is responsible to stop the adjoint tracking phase when:
//			-a)The adjoint track reaches the external surface.  
//			-b)The being tracked adjoint dynamic particle get an energy higher than the maximum energy of the external source.
//			-c)The adjoint track enters the volume delimited by the adjoint source.  
//		In the case a) the info (energy,weight,...) of the adjoint dynamic particle associated to the track
//			 when crossing the external source is registered and in the next event a forward primary is generated. In the other cases b) and c)
//			The next generated fwd particle is killed before being tracked and the next tracking of an adjoint particle is started directly. 	
//

#ifndef G4AdjointSteppingAction_h
#define G4AdjointSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4AdjointCrossSurfChecker;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4AdjointSteppingAction : public G4UserSteppingAction
{
  public:
    G4AdjointSteppingAction();
   ~G4AdjointSteppingAction();

    void UserSteppingAction(const G4Step*);
    
    inline void SetExtSourceEMax(G4double Emax){ext_sourceEMax=Emax;} 
    inline void SetStartEvent(G4bool aBool){start_event =aBool;}
    inline G4bool GetDidAdjParticleReachTheExtSource(){return did_adj_part_reach_ext_source;}
    inline G4ThreeVector GetLastMomentum(){return last_momentum;}
    inline G4ThreeVector GetLastPosition(){return last_pos;}
    inline G4double GetLastEkin(){return last_ekin;}
    inline G4double GetLastWeight(){return last_weight;}
    inline void SetPrimWeight(G4double weight){prim_weight=weight;} 
    inline G4ParticleDefinition* GetLastPartDef(){return last_part_def;}
    inline void SetUserAdjointSteppingAction( G4UserSteppingAction* anAction) {
    	                               theUserAdjointSteppingAction = anAction;}
    inline void SetUserForwardSteppingAction( G4UserSteppingAction* anAction) {
    	                                  theUserFwdSteppingAction = anAction;}
    inline void SetAdjointTrackingMode(G4bool aBool){is_adjoint_tracking_mode =aBool;}
    inline void ResetDidOneAdjPartReachExtSourceDuringEvent()
    						{did_one_adj_part_reach_ext_source_during_event =false;}
    inline void SetAdjointGeantinoTrackingMode(G4bool aBool){is_adjoint_geantino_tracking_mode =aBool;}
  private:

    G4double ext_sourceEMax;
    G4AdjointCrossSurfChecker* theG4AdjointCrossSurfChecker;
    G4bool start_event;
    
    G4bool did_adj_part_reach_ext_source;
    G4bool did_one_adj_part_reach_ext_source_during_event;
    G4ThreeVector last_momentum, last_pos; 
    G4double last_ekin;
    G4double last_weight ;
    G4double prim_weight ;
    G4ParticleDefinition* last_part_def;
    G4UserSteppingAction* theUserAdjointSteppingAction;
    G4UserSteppingAction* theUserFwdSteppingAction;
    G4bool is_adjoint_tracking_mode;
    G4bool is_adjoint_geantino_tracking_mode;


};
#endif

