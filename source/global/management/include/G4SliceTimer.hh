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
// G4SliceTimer
//
// Class description:
//
// Class for timer objects, able to measure elasped user/system process
// accumulated time.
//
// Note: Uses <sys/times.h> & <unistd.h> - POSIX.1 defined
//       If used, this header must be included in the source (.cc) file
//       and it must be the first header file to be included!

// Author: M.Asai, 23.10.06 - Derived from G4Timer implementation
// --------------------------------------------------------------------
#ifndef G4SLICE_TIMER_HH
#define G4SLICE_TIMER_HH 1

#if !(defined(WIN32) || defined(__MINGW32__))
#  include <sys/times.h>
#  include <unistd.h>
#else
#  include <time.h>
#  define _SC_CLK_TCK 1

extern "C"
{
  int sysconf(int);
};

// Structure returned by times()
//
struct tms
{
  clock_t tms_utime;  /* user time */
  clock_t tms_stime;  /* system time */
  clock_t tms_cutime; /* user time, children */
  clock_t tms_cstime; /* system time, children */
};

extern "C"
{
  extern clock_t times(struct tms*);
};
#endif /* WIN32 */

#include "G4Types.hh"
#include "G4ios.hh"

class G4SliceTimer
{
 public:
  G4SliceTimer();
  // Construct a timer object

  inline void Start();
  // Start timing
  inline void Stop();
  // Stop timing
  inline void Clear();
  // Clear accumulated times
  inline G4bool IsValid() const;
  // Return true if have a valid time (ie start() and stop() called)
  G4double GetRealElapsed() const;
  // Return the elapsed real time between last calling start() and stop()
  G4double GetSystemElapsed() const;
  // Return the elapsed system time between last calling start() and stop()
  G4double GetUserElapsed() const;
  // Return the elapsed user time between last calling start() and stop()

 private:
  clock_t fStartRealTime, fEndRealTime;
  // Real times (arbitrary time 0)
  tms fStartTimes, fEndTimes;
  // Timing structures (see times(2)) for start and end times

  G4double fRealElapsed = 0.0, fSystemElapsed = 0.0, fUserElapsed = 0.0;

  G4bool fValidTimes = true;
  // True after start and stop have both been called more than once and
  // an equal number of times
};

std::ostream& operator<<(std::ostream& os, const G4SliceTimer& t);
// Print the elapsed real,system and usertimes on os. Prints **s for times
// if !IsValid

#include "G4SliceTimer.icc"

#endif
