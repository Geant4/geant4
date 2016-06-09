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
// $Id: G4ParallelImportanceProcess.hh,v 1.12 2006/06/29 21:10:12 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceProcess
//
// Class description:
//
// Used internally by importance sampling in a "parallel" geometry.
// This is a G4ParallelTransport that also does importance
// sampling in the "parallel" geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceProcess_hh
#define G4ParallelImportanceProcess_hh G4ParallelImportanceProcess_hh

#include "G4ParallelTransport.hh"
#include "G4VTrackTerminator.hh"

class G4SamplingPostStepAction;
class G4VImportanceSplitExaminer;
class G4Nsplit_Weight;

class G4ParallelImportanceProcess : public G4ParallelTransport,
                                    public G4VTrackTerminator
{

public:  // with description

  G4ParallelImportanceProcess(const G4VImportanceSplitExaminer &aImportanceSplitExaminer,
                              G4VPGeoDriver &pgeodriver, 
                              G4VParallelStepper &aStepper,
                              const G4VTrackTerminator *TrackTerminator,
                              const G4String &aName = "ParallelImportanceProcess");  
    // initialise G4ParallelTransport and members

  virtual ~G4ParallelImportanceProcess();

  virtual G4VParticleChange *PostStepDoIt(const G4Track&,
                                          const G4Step&);
    // do the "parallel transport" and importance sampling.

  virtual void KillTrack() const;
    // used in case no scoring process follows that does the killing

  virtual const G4String &GetName() const;

private:

  G4ParallelImportanceProcess(const G4ParallelImportanceProcess &);
  G4ParallelImportanceProcess &operator=(const G4ParallelImportanceProcess &);
  
  virtual void Error(const G4String &m);

private:

  G4ParticleChange *fParticleChange;
  const G4VImportanceSplitExaminer &fImportanceSplitExaminer;  
  G4SamplingPostStepAction *fPostStepAction;
};

#endif
