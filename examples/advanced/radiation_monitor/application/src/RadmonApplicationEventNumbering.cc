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
//
// File name:     RadmonApplicationEventNumbering.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventNumbering.cc,v 1.3 2006/06/29 16:08:35 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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
