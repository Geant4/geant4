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
// G4AdjointTrackingAction
//
// Class description:
//
// This class represents actions taken place at the start/end point
// of processing one track during an adjoint simulation

// Author: L. Desorgher, SpaceIT GmbH
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// --------------------------------------------------------------------
#ifndef G4AdjointTrackingAction_hh
#define G4AdjointTrackingAction_hh 1

#include "globals.hh"
#include "G4UserTrackingAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4AdjointSteppingAction;
class G4Track;
class G4ParticleDefinition;

class G4AdjointTrackingAction : public G4UserTrackingAction
{

  public:

    G4AdjointTrackingAction(G4AdjointSteppingAction* anAction);
    virtual ~G4AdjointTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
    void RegisterAtEndOfAdjointTrack();
    void ClearEndOfAdjointTrackInfoVectors();

    inline void SetUserForwardTrackingAction(G4UserTrackingAction* anAction)
      { theUserFwdTrackingAction = anAction; }
    inline G4ThreeVector GetPositionAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_pos_vec[i]; }
    inline G4ThreeVector GetDirectionAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_direction_vec[i]; }
    inline G4double GetEkinAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_ekin_vec[i]; }
    inline G4double GetEkinNucAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_ekin_nuc_vec[i]; }
    inline G4double GetWeightAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_weight_vec[i]; }
    inline G4double GetCosthAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_cos_th_vec[i]; }
    inline const G4String& GetFwdParticleNameAtEndOfLastAdjointTrack()
      { return last_fwd_part_name; }
    inline G4int GetFwdParticlePDGEncodingAtEndOfLastAdjointTrack(std::size_t i=0)
      { return last_fwd_part_PDGEncoding_vec[i]; }
    inline G4bool GetIsAdjointTrackingMode()
      { return is_adjoint_tracking_mode; }
    inline G4int GetLastFwdParticleIndex(std::size_t i=0)
      { return last_fwd_part_index_vec[i]; }
    inline std::size_t GetNbOfAdointTracksReachingTheExternalSurface()
      { return last_pos_vec.size(); }
    inline void SetListOfPrimaryFwdParticles(std::vector<G4ParticleDefinition*>* aListOfParticles)
      { pListOfPrimaryFwdParticles = aListOfParticles; }

  private:

    G4AdjointSteppingAction* theAdjointSteppingAction = nullptr;
    G4UserTrackingAction* theUserFwdTrackingAction = nullptr;
    G4bool is_adjoint_tracking_mode = false;

    // Adjoint particle information on the external surface
    // ----------------------------------------------------
    G4ThreeVector last_pos;
    G4ThreeVector last_direction;
    G4double last_ekin = 0.0, last_ekin_nuc = 0.0;
      // last_ekin_nuc=last_ekin/nuc, nuc is 1 if not a nucleus
    G4double last_cos_th = 0.0;
    G4String last_fwd_part_name;
    G4int last_fwd_part_PDGEncoding = 0;
    G4double last_weight = 0.0;
    G4int last_fwd_part_index = 0;
    std::vector<G4ParticleDefinition*>* pListOfPrimaryFwdParticles = nullptr;

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
