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
// G4Timer class implementation
//
// Author: P.Kent, 21.08.95 - First implementation
// Revision: G.Cosmo, 29.04.97 - Added timings for Windows
// --------------------------------------------------------------------

#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4Exception.hh"

#include <iomanip>

#if defined(IRIX6_2)
#  if defined(_XOPEN_SOURCE) && (_XOPEN_SOURCE_EXTENDED == 1)
#    define __vfork vfork
#  endif
#endif

#ifdef WIN32
#  include <sys/types.h>
#  include <windows.h>

// extract milliseconds time unit
G4int sysconf(G4int a)
{
  if(a == _SC_CLK_TCK)
    return 1000;
  else
    return 0;
}

static clock_t filetime2msec(FILETIME* t)
{
  return (clock_t)((((G4float) t->dwHighDateTime) * 429496.7296) +
                   (((G4float) t->dwLowDateTime) * .0001));
}

clock_t times(struct tms* t)
{
  FILETIME ct = { 0, 0 }, et = { 0, 0 }, st = { 0, 0 }, ut = { 0, 0 },
           rt = { 0, 0 };
  SYSTEMTIME realtime;

  GetSystemTime(&realtime);
  SystemTimeToFileTime(&realtime, &rt);  // get real time in 10^-9 sec
  if(t != 0)
  {
    GetProcessTimes(GetCurrentProcess(), &ct, &et, &st, &ut);
    // get process time in 10^-9 sec
    t->tms_utime = t->tms_cutime = filetime2msec(&ut);
    t->tms_stime = t->tms_cstime = filetime2msec(&st);
  }
  return filetime2msec(&rt);
}
#endif /* WIN32 */

// Print timer status on std::ostream
//
std::ostream& operator<<(std::ostream& os, const G4Timer& t)
{
  // so fixed doesn't propagate
  std::stringstream ss;
  ss << std::fixed;
  if(t.IsValid())
  {
    ss << "User=" << t.GetUserElapsed() << "s Real=" << t.GetRealElapsed()
       << "s Sys=" << t.GetSystemElapsed() << "s";
#ifdef G4MULTITHREADED
    // avoid possible FPE error
    if(t.GetRealElapsed() > 1.0e-6)
    {
      G4double cpu_util = (t.GetUserElapsed() + t.GetSystemElapsed()) /
                          t.GetRealElapsed() * 100.0;
      ss << std::setprecision(1);
      ss << " [Cpu=" << std::setprecision(1) << cpu_util << "%]";
    }
#endif
  }
  else
  {
    ss << "User=****s Real=****s Sys=****s";
  }
  os << ss.str();

  return os;
}

G4double G4Timer::GetRealElapsed() const
{
  if(!fValidTimes)
  {
    G4Exception("G4Timer::GetRealElapsed()", "InvalidCondition", FatalException,
                "Timer not stopped or times not recorded!");
  }
  std::chrono::duration<G4double> diff = fEndRealTime - fStartRealTime;
  return diff.count();
}

G4double G4Timer::GetSystemElapsed() const
{
  if(!fValidTimes)
  {
    G4Exception("G4Timer::GetSystemElapsed()", "InvalidCondition",
                FatalException, "Timer not stopped or times not recorded!");
  }
  G4double diff = fEndTimes.tms_stime - fStartTimes.tms_stime;
  return diff / sysconf(_SC_CLK_TCK);
}

G4double G4Timer::GetUserElapsed() const
{
  if(!fValidTimes)
  {
    G4Exception("G4Timer::GetUserElapsed()", "InvalidCondition", FatalException,
                "Timer not stopped or times not recorded");
  }
  G4double diff = fEndTimes.tms_utime - fStartTimes.tms_utime;
  return diff / sysconf(_SC_CLK_TCK);
}
