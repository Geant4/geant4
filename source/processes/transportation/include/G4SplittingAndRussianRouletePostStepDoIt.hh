//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SplittingAndRussianRouletePostStepDoIt.hh,v 1.1 2003-08-19 16:37:23 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SplittingAndRussianRouletePostStepDoIt
//
// Class description:
//
// Used internally by importance and weight window sampling.
// Creates cloned tracks or kills tracks.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4SplittingAndRussianRouletePostStepDoIt_hh
#define G4SplittingAndRussianRouletePostStepDoIt_hh G4SplittingAndRussianRouletePostStepDoIt_hh

class G4VImportanceSplitExaminer;
class G4ParticleChange;
class G4Track;
class G4Step;
class G4Nsplit_Weight;
class G4VTrackTerminator;


class G4SplittingAndRussianRouletePostStepDoIt
{

public:  // with description

  explicit G4SplittingAndRussianRouletePostStepDoIt(const G4VTrackTerminator &
				    TrackTerminator);
    // simply construct

  ~G4SplittingAndRussianRouletePostStepDoIt();
    // simple destruct
  
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
