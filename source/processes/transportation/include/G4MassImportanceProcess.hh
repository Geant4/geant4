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
// $Id: G4MassImportanceProcess.hh,v 1.6 2002-10-16 16:26:58 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceProcess
//
// Class description:
//
// Used internally by importance sampling in the "mass" geometry.
// This process is a forced post step process. I will apply
// importance sampling if the particle crosses a boundary in the 
// "mass" geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceProcess_hh
#define G4MassImportanceProcess_hh G4MassImportanceProcess_hh

#include "G4VProcess.hh"
#include "G4ImportancePostStepDoIt.hh"
#include "G4VTrackTerminator.hh"
#include "G4ImportanceFinder.hh"

class G4VImportanceAlgorithm;
class G4VIStore;

class G4MassImportanceProcess : public G4VProcess, public G4VTrackTerminator
{

public:  // with description

  G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
			  const G4VIStore &aIstore,
			  const G4VTrackTerminator *TrackTerminator,
			  const G4String &aName = "MassImportanceProcess");
    // creates a G4ParticleChange

  virtual ~G4MassImportanceProcess();
    // delete the G4ParticleChange

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
    // make process beeing forced
  virtual G4VParticleChange *PostStepDoIt(const G4Track&, const G4Step&);
    // manage the importance sampling in the "mass" geometry

  virtual void KillTrack() const;
    // used in case no scoring process follows that does the killing

  virtual const G4String &GetName() const;


public:  // without description

  //  no operation in  AtRestDoIt and  AlongStepDoIt

  virtual G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*);
  virtual G4double 
  AtRestGetPhysicalInteractionLength(const G4Track& ,
				     G4ForceCondition*);
  
  virtual G4VParticleChange* 
  AtRestDoIt(const G4Track&, const G4Step&);


  virtual G4VParticleChange* 
  AlongStepDoIt(const G4Track&, const G4Step&);
  
private:
  
  G4MassImportanceProcess(const G4MassImportanceProcess &);
  G4MassImportanceProcess &operator=(const G4MassImportanceProcess &);
  
private:

  G4ParticleChange *fParticleChange;
  const G4VTrackTerminator *fTrackTerminator;
  const G4VImportanceAlgorithm &fImportanceAlgorithm;
  G4ImportanceFinder fImportanceFinder;
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
};

#endif
