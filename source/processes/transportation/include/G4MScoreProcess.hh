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
// $Id: G4MScoreProcess.hh,v 1.3 2002-04-10 13:14:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MScoreProcess
//
// Class description:
//
// Used internally by scoring in the "mass" world.
// This is a forced post step process messaging a "scorer" 
// derived from G4VPScorer.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MScoreProcess_hh 
#define G4MScoreProcess_hh G4MScoreProcess_hh 

#include "G4VProcess.hh"

class G4VPScorer;

class G4MScoreProcess : public G4VProcess
{

public:  // with description

  G4MScoreProcess(G4VPScorer &aScorer,
		  const G4String &aName = "MScoreProcess");
    // take reference to scorer and coppy particle name and
    // create a G4ParticleChange

  virtual ~G4MScoreProcess();
    // delete the G4ParticleChange

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
    // make processed being forced

  virtual G4VParticleChange * PostStepDoIt(const G4Track&, const G4Step&);
    // message "scorer" with  G4Step and a G4PStep from the "mass" 
    // geometry

public:  // without description

  // no operation in  AtRestDoIt and  AlongStepDoIt

  G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*) {return -1.0;}
  
  G4double 
  AtRestGetPhysicalInteractionLength(const G4Track&,
				     G4ForceCondition*) {return -1.0;}
  
  G4VParticleChange* AtRestDoIt(const G4Track&,
	                        const G4Step&) {return 0;}
  G4VParticleChange* AlongStepDoIt(const G4Track&,
				   const G4Step&) {return 0;}
  
private:

  G4MScoreProcess(const G4MScoreProcess &);
  G4MScoreProcess &operator=(const G4MScoreProcess &);

private:

  G4VPScorer &fScorer;  
};

#endif
