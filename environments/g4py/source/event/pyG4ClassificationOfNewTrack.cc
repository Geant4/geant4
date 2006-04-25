// $Id: pyG4ClassificationOfNewTrack.cc,v 1.2 2006-04-25 08:09:44 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ClassificationOfNewTrack.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ClassificationOfNewTrack.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ClassificationOfNewTrack()
{
 enum_<G4ClassificationOfNewTrack>("G4ClassificationOfNewTrack")
   .value("fUrgent",     fUrgent)
   .value("fWaiting",    fWaiting)
   .value("fPostpone",   fPostpone)
   .value("fKill",       fKill)
   .value("fWaiting_1",  fWaiting_1)
   .value("fWaiting_2",  fWaiting_2)
   .value("fWaiting_3",  fWaiting_3)
   .value("fWaiting_4",  fWaiting_4)
   .value("fWaiting_5",  fWaiting_5)
   .value("fWaiting_6",  fWaiting_6)
   .value("fWaiting_7",  fWaiting_7)
   .value("fWaiting_8",  fWaiting_8)
   .value("fWaiting_9",  fWaiting_9)
   .value("fWaiting_19", fWaiting_10)
   ;
}

