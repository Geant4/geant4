//
// File name:     RadmonRunAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunAction.cc,v 1.1 2005-11-24 02:31:56 capra Exp $
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
