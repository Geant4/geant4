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
/// \file PeriodicBoundaryProcess.hh
/// \brief Definition of the PeriodicBoundaryProcess class

/*
 * Based on 'G4pbc'.
 * Copyright (c) 2020 Amentum Pty Ltd
 * team@amentum.space
 * The original open-source version of this code
 * may be found at https://github.com/amentumspace/g4pbc
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * team@amentum.space
 */
#ifndef PeriodicBoundaryProcess_h
#define PeriodicBoundaryProcess_h 1

#include "ParticleChangeForPeriodic.hh"

#include "G4Electron.hh"
#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class G4Step;

enum ProcessStatus
{
  Undefined,
  Reflection,
  Cycling,
  StepTooSmall,
  NotAtBoundary
};

class PeriodicBoundaryProcess : public G4VDiscreteProcess
{
  public:
    explicit PeriodicBoundaryProcess(const G4String& processName = "PBC",
                                     G4ProcessType type = fNotDefined, bool per_x = true,
                                     bool per_y = true, bool per_z = true);

    ~PeriodicBoundaryProcess() override = default;

    PeriodicBoundaryProcess(const PeriodicBoundaryProcess& right) = delete;

    PeriodicBoundaryProcess& operator=(const PeriodicBoundaryProcess& right) = delete;

  public:
    G4bool IsApplicable(const G4ParticleDefinition&) override;

    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) override;

    ProcessStatus GetStatus() const;

    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  protected:
    ParticleChangeForPeriodic fParticleChange;

  private:
    void BoundaryProcessVerbose();
    std::map<ProcessStatus, G4String> fStatusMessages = {{Undefined, " *** Undefined *** "},
                                                         {NotAtBoundary, " *** NotAtBoundary *** "},
                                                         {Reflection, " *** Reflection *** "},
                                                         {Cycling, " *** Periodic *** "},
                                                         {StepTooSmall, " *** StepTooSmall *** "}};

    ProcessStatus fTheStatus = Undefined;
    G4ThreeVector fOldPosition;
    G4ThreeVector fNewPosition;
    G4ThreeVector fOldMomentum;
    G4ThreeVector fNewMomentum;
    G4ThreeVector fOldPolarization;
    G4ThreeVector fNewPolarization;
    G4ThreeVector fTheGlobalNormal;
    G4double fkCarTolerance;
    G4bool fPeriodicX = true;
    G4bool fPeriodicY = true;
    G4bool fPeriodicZ = true;
};

inline G4bool PeriodicBoundaryProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  G4bool applicable = false;
  // applied only for electrons
  if (&aParticleType == G4Electron::Electron()) {
    applicable = true;
  }
  return applicable;
}

inline ProcessStatus PeriodicBoundaryProcess::GetStatus() const
{
  return fTheStatus;
}

#endif
