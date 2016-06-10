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
// $Id: G4SamplingPostStepAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
// ----------------------------------------------------------------------
// Class G4SamplingPostStepAction
//
// Class description:
//
// Used internally by importance and weight window sampling.
// Creates cloned tracks or kills tracks.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4SamplingPostStepAction_hh
#define G4SamplingPostStepAction_hh G4SamplingPostStepAction_hh

class G4VImportanceSplitExaminer;
class G4ParticleChange;
class G4Track;
class G4Step;
class G4Nsplit_Weight;
class G4VTrackTerminator;

class G4SamplingPostStepAction
{

public:  // with description

  explicit G4SamplingPostStepAction(const G4VTrackTerminator &TrackTerminator);
    // Constructor

  ~G4SamplingPostStepAction();
    // Destructor
  
  void DoIt(const G4Track& aTrack, 
            G4ParticleChange *aParticleChange, 
            const G4Nsplit_Weight &nw);
    // Do the PostStepDoIt part common to importance and weight window
    // sampling in the 
    // "mass" and "parallel" geometry.
  
private:

  void Split(const G4Track &aTrack,
             const G4Nsplit_Weight &nw,
             G4ParticleChange *aParticleChange);

  const G4VTrackTerminator &fTrackTerminator;

};

#endif
