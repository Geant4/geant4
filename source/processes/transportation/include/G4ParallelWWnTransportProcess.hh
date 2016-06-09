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
// $Id: G4ParallelWWnTransportProcess.hh,v 1.5 2006/06/29 21:10:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4ParallelWWnTransportProcess
//
// Class description:
//
// Used internally by weight window technique in a "parallel" geometry.
// This is a G4ParallelTransport that also does weight window
// sampling in the "parallel" geometry. It is used in case
// the ww is applied only on the paralle boundaries.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelWWnTransportProcess_hh
#define G4ParallelWWnTransportProcess_hh G4ParallelWWnTransportProcess_hh

#include "G4ParallelTransport.hh"
#include "G4VTrackTerminator.hh"

class G4SamplingPostStepAction;
class G4VWeightWindowExaminer;
class G4Nsplit_Weight;

class G4ParallelWWnTransportProcess : public G4ParallelTransport,
                                      public G4VTrackTerminator
{

public:  // with description

  G4ParallelWWnTransportProcess(const G4VWeightWindowExaminer 
                                &aWeightWindowExaminer,
                                G4VPGeoDriver &pgeodriver, 
                                G4VParallelStepper &aStepper,
                                const G4VTrackTerminator *TrackTerminator,
                                const G4String &aName = 
                                "ParallelWWnTransportProcess");  
    // initialise G4ParallelTransport and members

  virtual ~G4ParallelWWnTransportProcess();


  virtual G4VParticleChange *PostStepDoIt(const G4Track&,
                                          const G4Step&);
    // do the "parallel transport" and weight window sampling

  virtual void KillTrack() const;
    // used in case no scoring process follows that does the killing

  virtual const G4String &GetName() const;


private:

  G4ParallelWWnTransportProcess(const G4ParallelWWnTransportProcess &);
  G4ParallelWWnTransportProcess &operator=(const G4ParallelWWnTransportProcess &);
  
  virtual void Error(const G4String &m);

private:

  G4ParticleChange *fParticleChange;
  const G4VWeightWindowExaminer &fWeightWindowExaminer;  
  G4SamplingPostStepAction *fPostStepAction;
};

#endif
