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
// $Id: $
//
//---------------------------------------------------------------
//
// G4BOptnCloning
//
// Class Description:
//    A G4VBiasingOperation to clone a track, allowing to set
//    weight arbitrary weights.
//    
//
//---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi


#ifndef G4BOptnCloning_hh
#define G4BOptnCloning_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"

class G4BOptnCloning : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnCloning(G4String name);
  // -- destructor:
  virtual ~G4BOptnCloning();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // -------------------------------------------
  // -- Unsed:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& ) {return 0;}
  virtual G4VParticleChange*                             ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
										 const G4Track*,
										 const G4Step*,
										 G4bool& ) {return 0;}
  // -- Used:
  virtual G4double                                  DistanceToApplyOperation( const G4Track*,
									      G4double,
									      G4ForceCondition* condition)
  {
    *condition = NotForced; return 0; // -- acts immediately
  }
  virtual G4VParticleChange*                    GenerateBiasingFinalState( const G4Track*,
									   const G4Step*  );
  
public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  void SetCloneWeights(G4double clone1Weight, G4double clone2Weight) {fClone1W = clone1Weight ; fClone2W = clone2Weight;}

  G4Track* GetCloneTrack() const { return fCloneTrack; }

private:
  G4double         fClone1W,
                   fClone2W;
  G4ParticleChange fParticleChange;
  G4Track*         fCloneTrack;
};

#endif
