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
// File name:     RadmonSteppingAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingAction.cc,v 1.2 2006-06-28 13:57:43 gunter Exp $
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
