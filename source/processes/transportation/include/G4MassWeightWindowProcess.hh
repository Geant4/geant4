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
// $Id: G4MassWeightWindowProcess.hh,v 1.4 2003/11/26 14:51:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// Class G4MassWeioghtWindowProcess
//
// Class description:
//
// Used internally by weight window technique in the "mass" geometry.
// This process is a forced post step process. It will apply
// weight window biasing on boundaries, collisions 
// or both according to the G4PlaceOfAction argument.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassWeioghtWindowProcess_hh
#define G4MassWeioghtWindowProcess_hh G4MassWeioghtWindowProcess_hh

#include "G4VProcess.hh"
#include "G4VTrackTerminator.hh"
#include "G4PlaceOfAction.hh"

class G4SamplingPostStepAction;
class G4VWeightWindowAlgorithm;
class G4VWeightWindowStore;

class G4MassWeightWindowProcess : public G4VProcess, public G4VTrackTerminator
{

public:  // with description

  G4MassWeightWindowProcess(const G4VWeightWindowAlgorithm &
                             aWeightWindowAlgorithm,
                             const G4VWeightWindowStore &aWWStore,
                             const G4VTrackTerminator *TrackTerminator,
                             G4PlaceOfAction placeOfAction,
                             const G4String &aName = 
                             "MassWeightWindowProcess");
    // creates a G4ParticleChange

  virtual ~G4MassWeightWindowProcess();
    // delete the G4ParticleChange

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                       G4double   previousStepSize,
                                       G4ForceCondition* condition);
    // make process beeing forced
  virtual G4VParticleChange *PostStepDoIt(const G4Track&, const G4Step&);
    // aply weight window sampling

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
  
  G4MassWeightWindowProcess(const G4MassWeightWindowProcess &);
  G4MassWeightWindowProcess &operator=(const G4MassWeightWindowProcess &);
  
private:

  G4ParticleChange *fParticleChange;
  const G4VWeightWindowAlgorithm &fWeightWindowAlgorithm;
  const G4VWeightWindowStore &fWeightWindowStore;
  G4SamplingPostStepAction *fPostStepAction;
  G4PlaceOfAction fPlaceOfAction;

};

#endif
