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
// $Id: G4StateManager.cc,v 1.8 2002-12-05 02:32:21 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      ---------------- G4StateManager ----------------
//             by Gabriele Cosmo, November 1996
// ------------------------------------------------------------

#include "G4StateManager.hh"

// Initialization of the static pointer of the single class instance
G4StateManager* G4StateManager::theStateManager = 0;

G4StateManager::G4StateManager()
 : theCurrentState(G4State_PreInit),
   thePreviousState(G4State_PreInit),
   theBottomDependent(0),
   suppressAbortion(0),
   msgptr(0)
{
}

G4StateManager::~G4StateManager()
{
  G4VStateDependent* state=0;

  while (theDependentsList.size()>0)
  {
    state = theDependentsList.back();
    theDependentsList.pop_back();
    for (G4std::vector<G4VStateDependent*>::iterator
         i=theDependentsList.begin(); i!=theDependentsList.end(); i++)
    {
      if (*i==state)
      {
	theDependentsList.erase(i);
	i--;
      }
    } 
    if ( state ) delete state;    
  } 
}

G4StateManager::G4StateManager(const G4StateManager &right)
  : theCurrentState(right.theCurrentState),
    thePreviousState(right.thePreviousState),
    theDependentsList(right.theDependentsList),
    theBottomDependent(right.theBottomDependent)
{
}

G4StateManager&
G4StateManager::operator=(const G4StateManager &right)
{
   if (&right == this) return *this;

   theCurrentState = right.theCurrentState;
   thePreviousState = right.thePreviousState;
   theDependentsList = right.theDependentsList;
   theBottomDependent = right.theBottomDependent;

   return *this;
}

G4int
G4StateManager::operator==(const G4StateManager &right) const
{
   return (this == &right);
}

G4int
G4StateManager::operator!=(const G4StateManager &right) const
{
   return (this != &right);
}


G4StateManager*
G4StateManager::GetStateManager()
{
    if (!theStateManager)
    {
      theStateManager = new G4StateManager;
    }
    return theStateManager;    
}

G4bool
G4StateManager::RegisterDependent(G4VStateDependent* aDependent, G4bool bottom)
{
    G4bool ack=true;
    if(!bottom)
    {
      theDependentsList.push_back(aDependent);
    }
    else
    { 
      if(theBottomDependent)
      {
        theDependentsList.push_back(theBottomDependent);
      }
      theBottomDependent = aDependent;
    }
    return ack;
}

G4bool
G4StateManager::DeregisterDependent(G4VStateDependent* aDependent)
{
  G4VStateDependent* tmp = 0;
  G4std::vector<G4VStateDependent*>::iterator i;
  for (i=theDependentsList.begin(); i!=theDependentsList.end(); i++)
    {
      if (**i==*aDependent) 
	{
	  tmp = *i;
	  theDependentsList.erase(i);
	} 
    }
  return (tmp != 0);
}

G4ApplicationState
G4StateManager::GetCurrentState() const
{
   return theCurrentState;
}

G4ApplicationState
G4StateManager::GetPreviousState() const
{
   return thePreviousState;
}

G4bool
G4StateManager::SetNewState(G4ApplicationState requestedState)
{ return SetNewState(requestedState,0); }

G4bool
G4StateManager::SetNewState(G4ApplicationState requestedState, const char* msg)
{
   if(requestedState==G4State_Abort && suppressAbortion>0) {
     if(suppressAbortion==2) return false;
     if(theCurrentState==G4State_EventProc) return false;
   }
   msgptr = msg;
   size_t i=0;
   G4bool ack = true;
   G4ApplicationState savedState = thePreviousState;
   thePreviousState = theCurrentState;
   theCurrentState = requestedState;
   while ((ack) && (i<theDependentsList.size()))
   {
     ack = theDependentsList[i]->Notify(requestedState);
     i++;
   }
   if(theBottomDependent)
   {
     ack = theBottomDependent->Notify(requestedState);
   }

   if(!ack)
   {
     theCurrentState = thePreviousState;
     thePreviousState = savedState;
   }
   msgptr = 0;
   return ack;
}

G4VStateDependent*
G4StateManager::RemoveDependent(const G4VStateDependent* aDependent)
{
  G4VStateDependent* tmp = 0;
  G4std::vector<G4VStateDependent*>::iterator i;
  for (i=theDependentsList.begin(); i!=theDependentsList.end(); i++)
    {
      if (**i==*aDependent) 
	{
	  tmp = *i;
	  theDependentsList.erase(i);
	} 
    }
  return tmp;
}

G4String
G4StateManager::GetStateString(G4ApplicationState aState) const
{
  G4String stateName;
  switch(aState)
  {
    case G4State_PreInit:
     stateName = "PreInit"; break;
    case G4State_Init:
     stateName = "Init"; break;
    case G4State_Idle:
     stateName = "Idle"; break;
    case G4State_GeomClosed:
     stateName = "GeomClosed"; break;
    case G4State_EventProc:
     stateName = "EventProc"; break;
    case G4State_Quit:
     stateName = "Quit"; break;
    case G4State_Abort:
     stateName = "Abort"; break;
    default:
     stateName = "Unknown";
  }
  return stateName;
}

//void G4StateManager::Pause()
//{
//  Pause("G4_pause> ");
//}
//
//void G4StateManager::Pause(const char* msg)
//{
//  G4String msgS = msg;
//  Pause(msgS);
//}
//
//void G4StateManager::Pause(G4String msg)
//{
//  G4UImanager::GetUIpointer()->PauseSession(msg);
//}
