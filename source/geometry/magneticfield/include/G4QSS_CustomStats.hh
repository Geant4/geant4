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
// QSSStats
//
// QSS statistics

// Authors:  Lucio Santi, Rodrigo Castro - 2018-2021
// --------------------------------------------------------------------
#ifndef _QSS_CUSTOM_STATS_HH_
#define _QSS_CUSTOM_STATS_HH_

#include <time.h>

#define GET_TIME(t0) (clock_gettime(CLOCK_MONOTONIC, &t0))
#define TIME_SECS(t0, t1) ((t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9)

// #define INTERPOLATE_ITERATIONS 1e2
// #define ESTIMATE_ITERATIONS 2
// #define ON_COMPUTE_STEP_ITERATIONS 20
// #define ON_COMPUTE_STEP_ITERATIONS_G4 2e3

#include "G4qss_misc.hh"
#include "G4Types.hh"

#include <atomic>

struct QSSStats
{
  G4double precision_dQMin;
  G4double precision_dQRel;
  G4int currentStep;
  G4int substeps;
  std::atomic<G4int> stepperSteps;
  std::map<G4int, std::map<G4int, G4int>> substepsByStepNumberByTrackID;

  G4double reset_time;
  G4double integration_time;

  G4int dqrel_changes[Qss_misc::VAR_IDX_END];
  G4int dqmin_changes[Qss_misc::VAR_IDX_END];
  G4double max_error[Qss_misc::VAR_IDX_END];

  QSSStats()
  {
    substeps = 0;
    reset_time = 0;
    integration_time = 0;

    for (size_t i = 0; i < Qss_misc::VAR_IDX_END; i++) {
      dqrel_changes[i] = 0;
      dqmin_changes[i] = 0;
      max_error[i] = 0;
    }
  };

  void print() const
  {
    G4int steps = stepperSteps.load();

    std::vector<std::string> vars{"x", "y", "z", "vx", "vy", "vz"};

    G4double avg_substeps = (G4double)substeps / steps;
    G4double avg_integration_time = (G4double)integration_time / steps;
    G4double avg_substeps_integration_time = (G4double)integration_time / substeps;
    G4double avg_reset_time = (G4double)reset_time / steps;

    std::stringstream ss;

    ss << "QSS stats:" << std::endl;
    ss << "dQMin: " << precision_dQMin << std::endl;
    ss << "dQRel: " << precision_dQRel << std::endl;

    ss << " Total steps: " << steps << std::endl
       << " Total substeps: " << substeps << std::endl
       << " Substeps average per step: " << avg_substeps << std::endl;

    ss << " Substeps by track-step:" << std::endl;
    for (auto it = substepsByStepNumberByTrackID.begin(); it != substepsByStepNumberByTrackID.end();
         ++it)
    {
      ss << "  Track #" << it->first << std::endl;
      for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        ss << "    Step " << it2->first << " => " << it2->second << " substeps" << std::endl;
      }
    }

    ss << " Integration time: " << integration_time << std::endl
       << " Integration time average (step): " << avg_integration_time << std::endl
       << " Integration time average (substep): " << avg_substeps_integration_time << std::endl;

    ss << " Reset time: " << reset_time << std::endl
       << " Reset time average: " << avg_reset_time << std::endl;

    for (G4int index = 0; index < Qss_misc::VAR_IDX_END; index++) {
      ss << " Variable " << vars[index] << ":" << std::endl;
      ss << "  dQRel changes: " << dqrel_changes[index] << std::endl;
      ss << "  dQMin changes: " << dqmin_changes[index] << std::endl;
      ss << "  Max error: " << max_error[index] << std::endl;
    }

    std::cout << ss.rdbuf();
  };
};

#endif
