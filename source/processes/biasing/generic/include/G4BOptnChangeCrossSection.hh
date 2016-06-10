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
// G4BOptnChangeCrossSection
//
// Class Description:
//    A G4VBiasingOperation to change a process cross-section.
//    
//
//---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi


#ifndef G4BOptnChangeCrossSection_hh
#define G4BOptnChangeCrossSection_hh 1

#include "G4VBiasingOperation.hh"
class G4InteractionLawPhysical;

class G4BOptnChangeCrossSection : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnChangeCrossSection(G4String name);
  // -- destructor:
  virtual ~G4BOptnChangeCrossSection();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // ----------------------------------------------
  // -- Used:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*,
										 G4ForceCondition& proposeForceCondition );
  
  // -- Unused:
  virtual G4VParticleChange*   ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
						       const G4Track*,
						       const G4Step*,
						       G4bool&          )  {return 0;}
  virtual G4double           DistanceToApplyOperation( const G4Track*,
						       G4double,
						       G4ForceCondition*)  {return DBL_MAX;}
  virtual G4VParticleChange* GenerateBiasingFinalState( const G4Track*,
						        const G4Step*   )  {return 0;}
  
public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  // -- return concrete type of interaction law:
  G4InteractionLawPhysical* GetBiasedExponentialLaw() {return fBiasedExponentialLaw;}
  // -- set biased cross-section:
  void     SetBiasedCrossSection(G4double xst);
  G4double GetBiasedCrossSection() const;
  // -- Sample underneath distribution:
  void Sample();
  // -- Update for a made step, without resampling:
  void UpdateForStep( G4double stepLength );
  // -- set/get if interaction occured. 
  // -- Interaction flag turned off in "Sample()" method.
  G4bool GetInteractionOccured() const { return fInteractionOccured; }
  void   SetInteractionOccured() { fInteractionOccured = true; }
  


private:
  G4InteractionLawPhysical* fBiasedExponentialLaw;
  G4bool                      fInteractionOccured;
};

#endif
