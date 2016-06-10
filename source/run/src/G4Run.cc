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
// $Id: G4Run.cc 70225 2013-05-27 10:10:15Z gcosmo $
//

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

G4Run::G4Run()
:runID(0),numberOfEvent(0),numberOfEventToBeProcessed(0),HCtable(0),DCtable(0)
{ eventVector = new std::vector<const G4Event*>; }

G4Run::~G4Run()
{
  // Objects made by local thread should not be deleted by the master thread
  G4RunManager::RMType rmType = G4RunManager::GetRunManager()->GetRunManagerType();
  if(rmType != G4RunManager::masterRM)
  {
    std::vector<const G4Event*>::iterator itr = eventVector->begin();
    for(;itr!=eventVector->end();itr++)
    { delete *itr; }
  }
  delete eventVector;
}

void G4Run::RecordEvent(const G4Event*)
{ numberOfEvent++; }

void G4Run::Merge(const G4Run* right)
{
  numberOfEvent += right->numberOfEvent; 
  std::vector<const G4Event*>::iterator itr = right->eventVector->begin();
  for(;itr!=right->eventVector->end();itr++)
  { eventVector->push_back(*itr); }
}

void G4Run::StoreEvent(G4Event* evt)
{ eventVector->push_back(evt); }

