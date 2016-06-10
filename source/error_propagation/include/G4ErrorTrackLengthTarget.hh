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
// $Id: G4ErrorTrackLengthTarget.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
// Class Description:
//
// G4ErrorTarget class: limits step when track length is bigger than
// theMaximumTrackLength. It is a G4VDiscreteProcess: limits the step in
// PostStepGetPhysicalInteractionLength().

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorTrackLengthTarget_hh
#define G4ErrorTrackLengthTarget_hh

#include "G4ios.hh" 
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Step.hh" 
#include "G4ErrorTarget.hh"

//---------------------------------------------------------------------------- 

class G4ErrorTrackLengthTarget : public G4VDiscreteProcess,
                                 public G4ErrorTarget
{
 public:  // with description

  G4ErrorTrackLengthTarget(const G4double maxTrkLength );
    // Constructs and add this process to G4ProcessManager of all particles
  virtual ~G4ErrorTrackLengthTarget(){}
  
    // These methods are dummy, as the step limitation is done in the
    // PostStepGetPhysicalInteractionLength().

  virtual G4double GetDistanceFromPoint( const G4ThreeVector&,
                                         const G4ThreeVector& ) const
    { return DBL_MAX; }   
  virtual G4double GetDistanceFromPoint( const G4ThreeVector& ) const
    { return DBL_MAX; } 

  virtual G4double
  PostStepGetPhysicalInteractionLength( const G4Track& track,
                                              G4double previousStepSize,
                                              G4ForceCondition* condition );
    // Checks if the maximum track length has been reached

  virtual  G4double GetMeanFreePath(const class G4Track & track,
                                    G4double, G4ForceCondition *);
    // Mean free path = theMaximumTrackLength - track.GetTrackLength()

   virtual G4VParticleChange* PostStepDoIt( const G4Track&, const G4Step& );
                       
   virtual void Dump( const G4String& msg ) const;

 private:

  G4double theMaximumTrackLength;
  G4VParticleChange theParticleChange;
};
   
#endif
