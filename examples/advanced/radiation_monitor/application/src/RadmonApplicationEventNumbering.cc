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
//
// File name:     RadmonApplicationEventNumbering.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventNumbering.cc,v 1.2 2006-06-28 13:45:48 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationEventNumbering.hh"
#include "G4Event.hh"

void                                            RadmonApplicationEventNumbering :: OnBeginOfEvent(const G4Event * event)
{
 if (dumpEvery==0)
  return;
  
 G4int eventId(event->GetEventID());

 if ((eventId%dumpEvery) == 0)
  G4cout << "RadmonApplicationEventNumbering::OnBeginOfEvent: Begin of event " << eventId << '.' << G4endl;
}
