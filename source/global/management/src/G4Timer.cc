// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Timer.cc,v 1.1 1999-01-07 16:09:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// class G4Timer
//
// Implementation
// 29.04.97 G.Cosmo Added timings for Windows/NT

#include "G4Timer.hh"
#include "G4ios.hh"

#ifdef WIN32
#  include <sys/types.h>
#  include <windows.h>

   // extract milliseconds time unit
   int sysconf(int a){
     if( a == _SC_CLK_TCK ) return 1000;
     else return 0;
   }

   static clock_t filetime2msec( FILETIME* t ){
   
      return (clock_t)((((float)t->dwHighDateTime)*429496.7296)+
              (((float)t->dwLowDateTime)*.0001) );
   }


   clock_t times(struct tms * t){
           FILETIME      ct = {0,0}, et = {0,0}, st = {0,0}, ut = {0,0}, rt = {0,0};
           SYSTEMTIME realtime;

           GetSystemTime( &realtime );
           SystemTimeToFileTime( &realtime, &rt ); // get real time in 10^-9 sec
           if( t != NULL ){
                   GetProcessTimes( GetCurrentProcess(), &ct, &et, &st, &ut);// get process time in 10^-9 sec
                   t->tms_utime = t->tms_cutime = filetime2msec(&ut);
                   t->tms_stime = t->tms_cstime = filetime2msec(&st);
           }
           return filetime2msec(&rt);
   }
#endif /* WIN32 */

// Print timer status n ostream
ostream& operator << (ostream& os, const G4Timer& t)
{
    if (t.IsValid())
	{
	    os << "User=" << t.GetUserElapsed()
	       << "s Real=" << t.GetRealElapsed()
	       << "s Sys=" << t.GetSystemElapsed() << "s";
	}
    else
	{
	    os << "User=****s Real=****s Sys=****s";
	}
    return os;
}

G4double G4Timer::GetRealElapsed() const
{
    if (!fValidTimes)
	{
	    G4Exception("G4Timer::GetRealElapsed - Timer not stopped or times not recorded");
	}
    G4double diff=fEndRealTime-fStartRealTime;
    return diff/sysconf(_SC_CLK_TCK);
}


G4double G4Timer::GetSystemElapsed() const
{
    if (!fValidTimes)
	{
	    G4Exception("G4Timer::GetSystemElapsed - Timer not stopped or times not recorded");
	}
    G4double diff=fEndTimes.tms_stime-fStartTimes.tms_stime;
    return diff/sysconf(_SC_CLK_TCK);
}

G4double G4Timer::GetUserElapsed() const
{
    if (!fValidTimes)
	{
	    G4Exception("G4Timer::GetUserElapsed - Timer not stopped or times not recorded");
	}
    G4double diff=fEndTimes.tms_utime-fStartTimes.tms_utime;
    return diff/sysconf(_SC_CLK_TCK);
}





