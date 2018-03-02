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
// $Id: G4StateManager.cc 108486 2018-02-15 14:47:25Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      ---------------- G4StateManager ----------------
//             by Gabriele Cosmo, November 1996
// ------------------------------------------------------------

#include "G4StateManager.hh"
#include "G4ios.hh"

// Initialization of the static pointer of the single class instance
//
G4ThreadLocal G4StateManager* G4StateManager::theStateManager = 0;
G4int G4StateManager::verboseLevel = 0;

G4StateManager::G4StateManager()
 : theCurrentState(G4State_PreInit),
   thePreviousState(G4State_PreInit),
   theBottomDependent(0),
   suppressAbortion(0),
   msgptr(0),
   exceptionHandler(0)
{
#ifdef G4MULTITHREADED
  G4iosInitialization();
#endif
}

G4StateManager::~G4StateManager()
{
  G4VStateDependent* state=0;

  while (theDependentsList.size()>0)
  {
    state = theDependentsList.back();
    theDependentsList.pop_back();
    for (std::vector<G4VStateDependent*>::iterator
         i=theDependentsList.begin(); i!=theDependentsList.end();)
    {
      if (*i==state)
      {
        i = theDependentsList.erase(i);
      }
      else
      {
        ++i;
      }
    } 
    if ( state )  { delete state; }
  } 
#ifdef G4MULTITHREADED_DEACTIVATE
  G4iosFinalization();
#endif
}

// -------------------------------------------------------------------------
// No matter how copy-constructor and operators below are implemented ...
// just dummy implementations, since not relevant for the singleton and
// declared private.
//
G4StateManager::G4StateManager(const G4StateManager &right)
  : theCurrentState(right.theCurrentState),
    thePreviousState(right.thePreviousState),
    theDependentsList(right.theDependentsList),
    theBottomDependent(right.theBottomDependent),
    suppressAbortion(right.suppressAbortion),
    msgptr(right.msgptr),
    exceptionHandler(right.exceptionHandler)
{
}

G4StateManager&
G4StateManager::operator=(const G4StateManager &right)
{
   if (&right == this)  { return *this; }

   theCurrentState = right.theCurrentState;
   thePreviousState = right.thePreviousState;
   theDependentsList = right.theDependentsList;
   theBottomDependent = right.theBottomDependent;
   suppressAbortion = right.suppressAbortion;
   msgptr = right.msgptr;
   exceptionHandler = right.exceptionHandler;

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
//
// -------------------------------------------------------------------------

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
  for (std::vector<G4VStateDependent*>::iterator i=theDependentsList.begin();
       i!=theDependentsList.end();)
  {
    if (**i==*aDependent) 
    {
      tmp = *i;
      i = theDependentsList.erase(i);
    }
    else
    {
      ++i;
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
   if(requestedState==G4State_Abort && suppressAbortion>0)
   {
     if(suppressAbortion==2)  { return false; }
     if(theCurrentState==G4State_EventProc)  { return false; }
   }
   msgptr = msg;
   size_t i=0;
   G4bool ack = true;
   G4ApplicationState savedState = thePreviousState;
   thePreviousState = theCurrentState;
   
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
   { thePreviousState = savedState; }
   else
   {
     theCurrentState = requestedState;
     if(verboseLevel>0)
     {
       G4cout<<"#### G4StateManager::SetNewState from "
             <<GetStateString(thePreviousState)<<" to "
             <<GetStateString(requestedState)<<G4endl;
     }
   }
   msgptr = 0;
   return ack;
}

G4VStateDependent*
G4StateManager::RemoveDependent(const G4VStateDependent* aDependent)
{
  G4VStateDependent* tmp = 0;
  for (std::vector<G4VStateDependent*>::iterator i=theDependentsList.begin();
       i!=theDependentsList.end();)
  {
    if (**i==*aDependent) 
    {
      tmp = *i;
      i = theDependentsList.erase(i);
    }
    else
    {
      ++i;
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
     stateName = "Unknown"; break;
  }
  return stateName;
}

void
G4StateManager::SetVerboseLevel(G4int val)
{
  verboseLevel = val;
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
