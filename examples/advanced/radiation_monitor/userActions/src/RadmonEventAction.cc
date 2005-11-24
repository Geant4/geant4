//
// File name:     RadmonEventAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonEventAction.cc,v 1.1 2005-11-24 02:31:56 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonEventAction.hh"
#include "RadmonEventActionObserver.hh"
#include "G4RunManager.hh"

void                                            RadmonEventAction :: AttachObserver(RadmonEventActionObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonEventAction :: DetachObserver(RadmonEventActionObserver * observer)
{
 observersSet.erase(observer);
}




void                                            RadmonEventAction :: BeginOfEventAction(const G4Event * event)
{
 ObserversSet::iterator i(observersSet.begin());
 const ObserversSet::iterator end(observersSet.end());
 
 while (i!=end)
 {
  (*i)->OnBeginOfEvent(event);
  i++;
 }
}



void                                            RadmonEventAction :: EndOfEventAction(const G4Event * event)
{
 ObserversSet::iterator i(observersSet.begin());
 const ObserversSet::iterator end(observersSet.end());
 
 while (i!=end)
 {
  (*i)->OnEndOfEvent(event);
  i++;
 }
}   





                                                RadmonEventAction :: RadmonEventAction()
{
 G4RunManager * runManager(G4RunManager::GetRunManager());
 
 runManager->SetUserAction(this);
}




RadmonEventAction *                             RadmonEventAction :: instance(0);
