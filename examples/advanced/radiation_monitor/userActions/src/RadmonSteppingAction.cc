//
// File name:     RadmonSteppingAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingAction.cc,v 1.1 2005-12-06 19:36:23 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonSteppingAction.hh"
#include "RadmonSteppingActionObserver.hh"
#include "G4RunManager.hh"

void                                            RadmonSteppingAction :: AttachObserver(RadmonSteppingActionObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonSteppingAction :: DetachObserver(RadmonSteppingActionObserver * observer)
{
 observersSet.erase(observer);
}




void                                            RadmonSteppingAction :: UserSteppingAction(const G4Step * step)
{
 ObserversSet::iterator i(observersSet.begin());
 const ObserversSet::iterator end(observersSet.end());
 
 while (i!=end)
 {
  (*i)->OnUserStepping(step);
  i++;
 }
}





                                                RadmonSteppingAction :: RadmonSteppingAction()
{
 G4RunManager * runManager(G4RunManager::GetRunManager());
 
 runManager->SetUserAction(this);
}




RadmonSteppingAction *                          RadmonSteppingAction :: instance(0);
