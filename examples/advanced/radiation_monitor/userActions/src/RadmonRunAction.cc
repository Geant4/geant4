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
// File name:     RadmonRunAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunAction.cc,v 1.2 2006-06-28 13:57:41 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonRunAction.hh"
#include "RadmonRunActionObserver.hh"
#include "G4RunManager.hh"

void                                            RadmonRunAction :: AttachObserver(RadmonRunActionObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonRunAction :: DetachObserver(RadmonRunActionObserver * observer)
{
 observersSet.erase(observer);
}




void                                            RadmonRunAction :: BeginOfRunAction(const G4Run * run)
{
 ObserversSet::iterator i(observersSet.begin());
 const ObserversSet::iterator end(observersSet.end());
 
 while (i!=end)
 {
  (*i)->OnBeginOfRun(run);
  i++;
 }
}



void                                            RadmonRunAction :: EndOfRunAction(const G4Run * run)
{
 ObserversSet::iterator i(observersSet.begin());
 const ObserversSet::iterator end(observersSet.end());
 
 while (i!=end)
 {
  (*i)->OnEndOfRun(run);
  i++;
 }
}   





                                                RadmonRunAction :: RadmonRunAction()
{
 G4RunManager * runManager(G4RunManager::GetRunManager());
 
 runManager->SetUserAction(this);
}




RadmonRunAction *                               RadmonRunAction :: instance(0);
