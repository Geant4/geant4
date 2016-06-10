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
// $Id: pyG4ClassificationOfNewTrack.cc 66892 2013-01-17 10:57:59Z gunter $
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

