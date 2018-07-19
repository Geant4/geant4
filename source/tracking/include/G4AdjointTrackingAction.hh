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
//
// $Id: G4AdjointTrackingAction.hh 75409 2013-11-01 00:31:00Z ldesorgh$
//
//
//---------------------------------------------------------------
//
// G4AdjointTrackingAction.hh
//
// class description:
//   This class represents actions taken place at the
//   the start/end point of processing one track during an adjoint simulation
//
//---------------------------------------------------------------



#ifndef G4AdjointTrackingAction_h
#define G4AdjointTrackingAction_h 1
#include "globals.hh"
#include "G4UserTrackingAction.hh"
#include "G4ThreeVector.hh"
#include <vector>
class G4AdjointSteppingAction;
class G4Track;
class G4ParticleDefinition;
///////////////////////////
class G4AdjointTrackingAction : public G4UserTrackingAction
///////////////////////////
{

//--------
public: // with description
//--------

// Constructor & Destructor
   G4AdjointTrackingAction(G4AdjointSteppingAction* anAction);
   virtual ~G4AdjointTrackingAction();

// Member functions
   virtual void PreUserTrackingAction(const G4Track*);
   virtual void PostUserTrackingAction(const G4Track*);
   void RegisterAtEndOfAdjointTrack();
   void ClearEndOfAdjointTrackInfoVectors();

//inline methods
   inline void SetUserForwardTrackingAction(G4UserTrackingAction* anAction){
	   	   	               theUserFwdTrackingAction = anAction;}
   inline G4ThreeVector GetPositionAtEndOfLastAdjointTrack(size_t i=0){ return last_pos_vec[i];}
   inline G4ThreeVector GetDirectionAtEndOfLastAdjointTrack(size_t i=0){ return last_direction_vec[i];}
   inline G4double GetEkinAtEndOfLastAdjointTrack(size_t i=0){ return last_ekin_vec[i];}
   inline G4double GetEkinNucAtEndOfLastAdjointTrack(size_t i=0){ return last_ekin_nuc_vec[i];}
   inline G4double GetWeightAtEndOfLastAdjointTrack(size_t i=0){return last_weight_vec[i];}
   inline G4double GetCosthAtEndOfLastAdjointTrack(size_t i=0){return last_cos_th_vec[i];}
   inline const G4String& GetFwdParticleNameAtEndOfLastAdjointTrack(){return last_fwd_part_name;}
   inline G4int GetFwdParticlePDGEncodingAtEndOfLastAdjointTrack(size_t i=0){return last_fwd_part_PDGEncoding_vec[i];}
   inline G4bool GetIsAdjointTrackingMode(){return is_adjoint_tracking_mode;}
   inline G4int GetLastFwdParticleIndex(size_t i=0){
             return last_fwd_part_index_vec[i];};
   inline size_t GetNbOfAdointTracksReachingTheExternalSurface(){return last_pos_vec.size();}
   inline void SetListOfPrimaryFwdParticles( std::vector<G4ParticleDefinition*>*
            aListOfParticles){pListOfPrimaryFwdParticles=aListOfParticles;}


private:
   G4AdjointSteppingAction* theAdjointSteppingAction;
   G4UserTrackingAction* theUserFwdTrackingAction;
   G4bool is_adjoint_tracking_mode;


   //adjoint particle information on the external surface
   //-----------------------------
   G4ThreeVector last_pos;
   G4ThreeVector last_direction;
   G4double last_ekin,last_ekin_nuc; //last_ekin_nuc=last_ekin/nuc, nuc is 1 if not a nucleus
   G4double last_cos_th;
   G4String last_fwd_part_name;
   G4int  last_fwd_part_PDGEncoding;
   G4double last_weight;
   G4int  last_fwd_part_index;
   std::vector<G4ParticleDefinition*>* pListOfPrimaryFwdParticles;

   std::vector<G4ThreeVector> last_pos_vec;
   std::vector<G4ThreeVector> last_direction_vec;
   std::vector<G4double>  last_ekin_vec;
   std::vector<G4double>  last_ekin_nuc_vec;
   std::vector<G4double>  last_cos_th_vec;
   std::vector<G4double> last_weight_vec;
   std::vector<G4int> last_fwd_part_PDGEncoding_vec;
   std::vector<G4int> last_fwd_part_index_vec;

};

#endif


