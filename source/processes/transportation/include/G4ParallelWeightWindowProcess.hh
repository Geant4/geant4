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
// $Id: G4ParallelWeightWindowProcess.hh,v 1.9 2003/11/26 14:51:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// ----------------------------------------------------------------------
// Class G4ParallelWeightWindowProcess
//
// Class description:
//
// Used internally by weight window technique in a "parallel" geometry.
// This process is a forced post step process. It will apply
// weight window biasing on collisions or on collisions and
// boundaries according to the G4PlaceOfAction argument. 
// It is not used in the case where it should apply only on boundaries.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelWeightWindowProcess_hh
#define G4ParallelWeightWindowProcess_hh G4ParallelWeightWindowProcess_hh

#include "G4VProcess.hh"
#include "G4VTrackTerminator.hh"
#include "G4PlaceOfAction.hh"

class G4SamplingPostStepAction;
class G4VWeightWindowExaminer;
class G4Nsplit_Weight;
class G4VParallelStepper;
class G4VParallelStepper;

class G4ParallelWeightWindowProcess : public G4VProcess,\
                                      public G4VTrackTerminator
{

public:  // with description

  G4ParallelWeightWindowProcess(const G4VWeightWindowExaminer 
                                &aWeightWindowExaminer,
                                G4VParallelStepper &aStepper,
                                const G4VTrackTerminator *TrackTerminator,
                                G4PlaceOfAction placeOfAction,
                                const G4String &aName = 
                                "ParallelWeightWindowProcess");  
    // initialise G4ParallelTransport and members

  virtual ~G4ParallelWeightWindowProcess();

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                       G4double   previousStepSize,
                                       G4ForceCondition* condition);
    // make process beeing forced
  virtual G4VParticleChange *PostStepDoIt(const G4Track&, const G4Step&);
    // apply weight window sampliing

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

  G4ParallelWeightWindowProcess(const G4ParallelWeightWindowProcess &);
  G4ParallelWeightWindowProcess &operator=(const G4ParallelWeightWindowProcess &);
  
  virtual void Error(const G4String &m);

private:

  G4ParticleChange *fParticleChange;
  const G4VWeightWindowExaminer &fWeightWindowExaminer;
  G4VParallelStepper &fStepper;
  G4SamplingPostStepAction *fPostStepAction;
  G4PlaceOfAction fPlaceOfAction;
};

#endif
