//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pyG4ClassificationOfNewTrack.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
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

