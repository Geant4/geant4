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
// G4Timer
//
// Class description:
//
// Class for timer objects, able to measure elasped user/system process time.
//
// Note: Uses <sys/times.h> & <unistd.h> - POSIX.1 defined
//       If used, this header must be included in the source (.cc) file and it
//       must be the first header file to be included!
//
// Member functions:
//
// G4Timer()
//   Construct a timer object
// Start()
//   Start timing
// Stop()
//   Stop timing
// G4bool IsValid()
//   Return true if have a valid time (ie start() and stop() called)
// G4double GetRealElapsed()
//   Return the elapsed real time between last calling start() and stop()
// G4double GetSystemElapsed()
//   Return the elapsed system time between last calling start() and stop()
// G4double GetUserElapsed()
//   Return the elapsed user time between last calling start() and stop()
//
// Operators:
//
// std::ostream& operator << (std::ostream& os, const G4Timer& t);
//   Print the elapsed real,system and usertimes on os. Prints **s for times
//   if !IsValid
//
// Member data:
//
// G4bool fValidTimes
//   True after start and stop have both been called more than once and
//   an equal number of times
// clock_t fStartRealTime,fEndRealTime
//   Real times (arbitrary time 0)
// tms fStartTimes,fEndTimes
//   Timing structures (see times(2)) for start and end times

// Author: P.Kent, 21.08.95 - First implementation
// Revision: G.Cosmo, 29.04.97 - Added timings for Windows
// --------------------------------------------------------------------
#ifndef G4TIMER_HH
#define G4TIMER_HH 1

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

#include <chrono>

class G4Timer
{
 public:
  inline void Start();
  inline void Stop();
  inline G4bool IsValid() const;
  inline const char* GetClockTime() const;
  G4double GetRealElapsed() const;
  G4double GetSystemElapsed() const;
  G4double GetUserElapsed() const;

 private:
  G4bool fValidTimes{false};
  using clock_type = std::chrono::high_resolution_clock;
  std::chrono::time_point<clock_type> fStartRealTime, fEndRealTime;
  tms fStartTimes, fEndTimes;
};

std::ostream& operator<<(std::ostream& os, const G4Timer& t);

#include "G4Timer.icc"

#endif
