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
// $Id: G4MultiEventAction.cc 90212 2016-01-27 18:33:12Z adotti $
//
//---------------------------------------------------------------
//
// G4MultiEventAction.cc
//
//   Created on: Jan 17, 2016
//       Author: adotti
//
// ---------------------------------------------------------------
//

#include "G4MultiEventAction.hh"
#include <algorithm>

void G4MultiEventAction::SetEventManager(G4EventManager* mgr)
{
  std::for_each( begin() , end() ,
      [mgr](G4UserEventActionUPtr& e) { e->SetEventManager(mgr); }
  );
}

void G4MultiEventAction::BeginOfEventAction(const G4Event* evt)
{
  std::for_each( begin() , end() ,
      [evt](G4UserEventActionUPtr& e) { e->BeginOfEventAction(evt); }
  );
}

void G4MultiEventAction::EndOfEventAction(const G4Event* evt)
{
  std::for_each( begin() , end() ,
      [evt](G4UserEventActionUPtr& e) { e->EndOfEventAction(evt); }
  );
}

