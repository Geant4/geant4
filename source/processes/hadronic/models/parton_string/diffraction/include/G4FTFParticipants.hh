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
//

#ifndef G4FTFParticipants_h
#define G4FTFParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
// ------------------------------------------------------------

#include "G4VParticipants.hh"
#include "G4FTFParameters.hh"
#include <vector>
#include "G4Nucleon.hh"
#include "G4V3DNucleus.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4InteractionContent.hh"


class G4FTFParticipants : public G4VParticipants {
  public:
    G4FTFParticipants();
    ~G4FTFParticipants();

    const G4FTFParticipants& operator=( const G4FTFParticipants& right ) = delete;
    G4bool operator==( const G4FTFParticipants& right ) const = delete;
    G4bool operator!=( const G4FTFParticipants& right ) const = delete;
    G4FTFParticipants( const G4FTFParticipants& right ) = delete;

    // New methods to get/set the impact parameter.
    // (Note: to get the impact parameter in fermi units, do:
    //          GetImpactParameter() / fermi )
    // By default, SampleBinInterval() is false and therefore the sampling of
    // the impact parameter is done as usual by FTF taking into account the
    // radius of target nucleus, and, in the case of nucleus-nucleus collisions,
    // also of radius of the projectile nucleus; if, instead, SampleBinInterval()
    // is true, then the square of the impact parameter is drawn from a flat
    // distribution in the specified interval [Bmin2, Bmax2].
    void SetImpactParameter( const G4double b_value );
    G4double GetImpactParameter() const;
    void SetBminBmax( const G4double bmin_value, const G4double bmax_value );
    G4bool SampleBinInterval() const;
    G4double GetBmin2() const;  // Minimum value of the square of the impact parameter
    G4double GetBmax2() const;  // Maximum value of the square of the impact parameter

    void GetList( const G4ReactionProduct& thePrimary, G4FTFParameters* theParameters );
    void StartLoop();
    G4bool Next();
    void SortInteractionsIncT();
    void ShiftInteractionTime();
    G4InteractionContent& GetInteraction();  
    void Clean();

  private:
    G4double Bimpact;
    G4bool   BinInterval;
    G4double Bmin2;
    G4double Bmax2;

    std::vector< G4InteractionContent* > theInteractions;
    G4int currentInteraction;
};


inline void G4FTFParticipants::SetImpactParameter( const G4double b_value ) {
  Bimpact = b_value;
}

inline G4double G4FTFParticipants::GetImpactParameter() const {
  return Bimpact;
}

inline void G4FTFParticipants::SetBminBmax( const G4double bmin_value, const G4double bmax_value ) {
  BinInterval = false;
  if ( bmin_value < 0.0 || bmax_value < 0.0 || bmax_value < bmin_value ) return;
  BinInterval = true;
  Bmin2 = bmin_value * bmin_value;
  Bmax2 = bmax_value * bmax_value;
}

inline G4bool G4FTFParticipants::SampleBinInterval() const {
  return BinInterval;
}

inline G4double G4FTFParticipants::GetBmin2() const {
  return Bmin2;
}

inline G4double G4FTFParticipants::GetBmax2() const {
  return Bmax2;
}

inline void G4FTFParticipants::StartLoop() {
  currentInteraction = -1;
}

inline G4bool G4FTFParticipants::Next() {
  return ++currentInteraction < static_cast< G4int >( theInteractions.size() );
}

inline G4InteractionContent& G4FTFParticipants::GetInteraction() {
  return *theInteractions[ currentInteraction ];
}

#endif

