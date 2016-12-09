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
// $Id: G4AdjointForcedInteractionForGamma.hh 80314 2014-04-10 12:23:52Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointForcedInteractionForGamma
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	12 September 2016 creation by L. Desorgher
//
//-------------------------------------------------------------
//	Documentation:
//		Class for the forced interaction of reverse gamma
//

#ifndef G4AdjointForcedInteractionForGamma_h
#define G4AdjointForcedInteractionForGamma_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include"G4PhysicsOrderedFreeVector.hh"


class G4PhysicsTable;
class G4Region;
class G4VParticleChange;
class G4ParticleChange;
class G4Track;
class G4VEmAdjointModel;
class G4AdjointCSMatrix;
class G4AdjointCSManager;
class G4Material;
class G4MaterialCutsCouple;
class G4Navigator;
class G4AdjointForcedInteractionForGamma : public G4VContinuousDiscreteProcess
{

public:

  G4AdjointForcedInteractionForGamma(G4String process_name);

  virtual ~G4AdjointForcedInteractionForGamma();
  
public:
  void PreparePhysicsTable(const G4ParticleDefinition&);
  void BuildPhysicsTable(const G4ParticleDefinition&);
  virtual G4double PostStepGetPhysicalInteractionLength(
                              const G4Track& track,
 			     G4double   previousStepSize,
 			     G4ForceCondition* condition);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
  virtual G4VParticleChange* AlongStepDoIt(const G4Track& track,const G4Step& step);
  inline void RegisterAdjointComptonModel(G4VEmAdjointModel* aAdjointComptonModel){theAdjointComptonModel = aAdjointComptonModel;}
  inline void RegisterAdjointBremModel(G4VEmAdjointModel* aAdjointBremModel){theAdjointBremModel = aAdjointBremModel;}
 
protected :// with description  

   virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);
   virtual G4double GetContinuousStepLimit(const G4Track& aTrack,
                                G4double  previousStepSize,
                                G4double  currentMinimumStep,
   			     G4double& currentSafety
                                                                );

private:
   G4VEmAdjointModel* theAdjointComptonModel;
   G4VEmAdjointModel* theAdjointBremModel;

   G4ParticleChange* fParticleChange;
   G4AdjointCSManager* theAdjointCSManager;
   
private:
  G4double lastAdjCS,lastFwdCS;

  G4int trackid;
  G4int nstep;  
  G4bool is_free_flight_gamma;
  G4bool copy_gamma_for_forced_interaction;
  G4int last_free_flight_trackid;

  G4double acc_track_length;
  G4double total_acc_nb_adj_interaction_length;
  G4double total_acc_nb_fwd_interaction_length;


  G4double acc_nb_adj_interaction_length;
  G4double acc_nb_fwd_interaction_length;
  G4bool continue_gamma_as_new_free_flight;
};  

#endif

