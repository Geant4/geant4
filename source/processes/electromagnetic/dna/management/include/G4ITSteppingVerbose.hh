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
  ~G4ITSteppingVerbose() override;

  // methods to be invoked in the SteppingManager
  void NewStep() override;
  void StepInfoForLeadingTrack() override;

  void AtRestDoItInvoked() override;
  void AtRestDoItOneByOne() override;

  void AlongStepDoItAllDone() override;
  void AlongStepDoItOneByOne() override;

  void PostStepDoItAllDone() override;
  void PostStepDoItOneByOne() override;

  void StepInfo() override;
  void TrackingStarted(G4Track*) override;
  void TrackingEnded(G4Track*) override;

  void DoItStarted() override;
  void PreStepVerbose(G4Track* track) override;
  void PostStepVerbose(G4Track* track) override;

  void DPSLStarted() override;
  void DPSLUserLimit() override;
  void DPSLPostStep() override;
  void DPSLAlongStep() override;
//  void DPSLAlongStepDoItOneByOne();
//  void DPSLPostStepDoItOneByOne();
  void VerboseTrack() override;
  void VerboseParticleChange() override;
  void ShowStep() const;
//

};

#endif

