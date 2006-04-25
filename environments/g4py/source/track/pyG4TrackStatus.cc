// $Id: pyG4TrackStatus.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4TrackStatus.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TrackStatus.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4TrackStatus()
{
  enum_<G4TrackStatus>("G4TrackStatus")
    .value("fAlive",                   fAlive)
    .value("fStopButAlive",            fStopButAlive)
    .value("fStopAndKill",             fStopAndKill)
    .value("fKillTrackAndSecondaries", fKillTrackAndSecondaries)
    .value("fSuspend",                 fSuspend)
    .value("fPostponeToNextEvent",     fPostponeToNextEvent)
    ;
}
