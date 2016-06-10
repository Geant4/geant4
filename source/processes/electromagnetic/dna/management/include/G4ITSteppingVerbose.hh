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
// $Id$
//
//---------------------------------------------------------------
//
// G4SteppingVerbose.hh
//
// class description:
//   This class manages the verbose outputs in G4ITStepProcessor.
//
//---------------------------------------------------------------

class G4SteppingVerbose;

#ifndef G4SteppingVerose_h
#define G4SteppingVerose_h 1

#include "G4VITSteppingVerbose.hh"

class G4ITSteppingVerbose : public G4VITSteppingVerbose
{
public:
  G4ITSteppingVerbose();
  ~G4ITSteppingVerbose();

  // methods to be invoked in the SteppingManager
  void NewStep();
  void StepInfoForLeadingTrack();

  void AtRestDoItInvoked();
  void AtRestDoItOneByOne();

  void AlongStepDoItAllDone();
  void AlongStepDoItOneByOne();

  void PostStepDoItAllDone();
  void PostStepDoItOneByOne();

  void StepInfo();
  void TrackingStarted(G4Track*);
  void TrackingEnded(G4Track*);

  void DoItStarted();
  void PreStepVerbose(G4Track* track);
  void PostStepVerbose(G4Track* track);

  void DPSLStarted();
  void DPSLUserLimit();
  void DPSLPostStep();
  void DPSLAlongStep();
//  void DPSLAlongStepDoItOneByOne();
//  void DPSLPostStepDoItOneByOne();
  void VerboseTrack();
  void VerboseParticleChange();
  void ShowStep() const;
//

};

#endif

