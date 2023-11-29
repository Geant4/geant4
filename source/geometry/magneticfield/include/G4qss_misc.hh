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
// Authors:  Lucio Santi, Rodrigo Castro         2018-21 

#ifndef _QSS_MISC_H_
#define _QSS_MISC_H_

typedef struct QSS_simulator_ *QSS_simulator;
typedef struct QSSSubstep_ *QSSSubstep;

namespace Qss_misc { 
  // Convention of Geant4 notation of indices
  constexpr unsigned int PXidx= 0;
  constexpr unsigned int PYidx= 1;
  constexpr unsigned int PZidx= 2;

  constexpr unsigned int VXidx= 3;
  constexpr unsigned int VYidx= 4;
  constexpr unsigned int VZidx= 5;

  // Method parameters & constants
  constexpr unsigned int MAX_QSS_STEPPER_ORDER= 3;
  constexpr unsigned int VAR_IDX_END= 6;
  constexpr unsigned int MIN_SUBSTEPS= 20;

  constexpr G4double INF= 1.0e20;
}

#if defined(WIN32) || defined(__MINGW32__)
#define unlikely(x)   (x)  // Until C++20 can be assumed
#define   likely(x)   (x)  //    >> ditto >>
#else
#define unlikely(x)   __builtin_expect((x),0)   // gcc/clang extension - not portable
#define   likely(x)   __builtin_expect((x),1)
#endif

// #define likely(x)   (x)  // [[likely]]     // The C++20 portable way
// #define likely(x)   (x)  // [[unlikely]]   //     >>    >> 
// This syntax appears to be part of C++20
// See
// - https://en.cppreference.com/w/cpp/language/attributes/likely
// - https://stackoverflow.com/questions/51797959/how-to-use-c20s-likely-unlikely-attribute-in-if-else-statement
// - https://usingstdcpp.org/2018/03/18/jacksonville18-iso-cpp-report/

#define SUBSTEP_STRUCT(sim, i)     (sim->substeps[i])
#define SUBSTEP_START(sim, i)      (sim->substeps[(i)].start_time)
#define SUBSTEP_X(sim, i)          (sim->substeps[(i)].x)
#define SUBSTEP_TX(sim, i)         (sim->substeps[(i)].tx)
#define SUBSTEP_LEN(sim, i)        (sim->substeps[(i)].len)

#define LAST_SUBSTEP_STRUCT(sim)     (SUBSTEP_STRUCT(sim, sim->cur_substep_idx-1))

#define CUR_SUBSTEP_START(sim)      (SUBSTEP_START(sim, sim->cur_substep_idx))
#define CUR_SUBSTEP_X(sim)          (SUBSTEP_X(sim, sim->cur_substep_idx))
#define CUR_SUBSTEP_TX(sim)         (SUBSTEP_TX(sim, sim->cur_substep_idx))
#define CUR_SUBSTEP_LEN(sim)        (SUBSTEP_LEN(sim, sim->cur_substep_idx))

#define CUR_SUBSTEP(sim)            (sim->cur_substep_idx)
#define LAST_SUBSTEP(sim)           (sim->cur_substep_idx-1)
#define MAX_SUBSTEP(sim)            (sim->max_substep_idx)
#define SUBSTEPS(sim)               (sim->substeps)

struct QSSSubstep_
{
  double x[Qss_misc::VAR_IDX_END*(Qss_misc::MAX_QSS_STEPPER_ORDER+1)];
  double tx[Qss_misc::VAR_IDX_END];

  double start_time;
  double len;
};

struct QSS_simulator_
{
  double x[Qss_misc::VAR_IDX_END*(Qss_misc::MAX_QSS_STEPPER_ORDER+1)];
  double tx[Qss_misc::VAR_IDX_END];

  double q[Qss_misc::VAR_IDX_END*(Qss_misc::MAX_QSS_STEPPER_ORDER+1)];
  double tq[Qss_misc::VAR_IDX_END];

  double nextStateTime[Qss_misc::VAR_IDX_END];
  double time;
  int minIndex;

  double dQMin[Qss_misc::VAR_IDX_END];
  double dQRel[Qss_misc::VAR_IDX_END];
  double lqu[Qss_misc::VAR_IDX_END];

  double alg[Qss_misc::VAR_IDX_END];
  double it;

  int *SD[Qss_misc::VAR_IDX_END];
  int states;

  QSSSubstep substeps;
  int cur_substep_idx;
  int max_substep_idx;
};

#endif
