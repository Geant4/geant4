// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StateManager.cc,v 1.2 1999-12-15 14:50:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4StateManager ----------------
//             by Gabriele Cosmo, November 1996
// ------------------------------------------------------------

#include "G4StateManager.hh"
#include "G4UImanager.hh"

// Initialization of the static pointer of the single class instance
G4StateManager* G4StateManager::theStateManager = 0;

G4StateManager::G4StateManager()
:theCurrentState(PreInit),thePreviousState(PreInit) {;}

G4StateManager::~G4StateManager()
{
   theDependentsList.clearAndDestroy();
}

G4StateManager::G4StateManager(const G4StateManager &right)
{
   *this = right;
}

G4StateManager& G4StateManager::operator=(const G4StateManager &right)
{
   *this = right;
   return *this;
}

G4int G4StateManager::operator==(const G4StateManager &right) const
{
   return (this == (G4StateManager *) &right);
}

G4int G4StateManager::operator!=(const G4StateManager &right) const
{
   return (this != (G4StateManager *) &right);
}


G4StateManager* G4StateManager::GetStateManager()
{
    if (!theStateManager)
        {
            theStateManager = new G4StateManager;
        }
    return theStateManager;    
}

G4bool G4StateManager::RegisterDependent(G4VStateDependent* aDependent)
{
   G4bool ack=true;
   theDependentsList.insert(aDependent);
   return ack;
}

G4bool G4StateManager::DeregisterDependent(G4VStateDependent* aDependent)
{
   G4VStateDependent* aD = theDependentsList.remove(aDependent);
   return (aD != NULL);
}

G4ApplicationState G4StateManager::GetCurrentState()
{
   return theCurrentState;
}

G4ApplicationState G4StateManager::GetPreviousState()
{
   return thePreviousState;
}

G4bool G4StateManager::SetNewState(G4ApplicationState requestedState)
{
   G4int i=0;
   G4bool ack = true;
   G4ApplicationState savedState = thePreviousState;
   thePreviousState = theCurrentState;
   theCurrentState = requestedState;
   while ((ack) && (i<theDependentsList.entries())) {
     ack = theDependentsList(i)->Notify(requestedState);
     i++;
   }

   if(!ack)
   {
     theCurrentState = thePreviousState;
     thePreviousState = savedState;
   }
   return ack;
}

G4VStateDependent* G4StateManager::RemoveDependent(const G4VStateDependent* aDependent)
{
   return theDependentsList.remove(aDependent);
}

G4String G4StateManager::GetStateString(G4ApplicationState aState)
{
  G4String stateName;
  switch(aState)
  {
    case PreInit:
     stateName = "PreInit"; break;
    case Init:
     stateName = "Init"; break;
    case Idle:
     stateName = "Idle"; break;
    case GeomClosed:
     stateName = "GeomClosed"; break;
    case EventProc:
     stateName = "EventProc"; break;
    case Quit:
     stateName = "Quit"; break;
    default:
     stateName = "Unknown";
  }
  return stateName;
}

void G4StateManager::Pause()
{
  Pause("G4_pause> ");
}

void G4StateManager::Pause(char* msg)
{
  G4String msgS = msg;
  Pause(msgS);
}

void G4StateManager::Pause(G4String msg)
{
  G4UImanager::GetUIpointer()->PauseSession(msg);
//  if(theCurrentState==EventProc)
//  {
//    G4ApplicationState savedState = thePreviousState;
//    SetNewState(Pause);
//    theCurrentState = thePreviousState;
//    thePreviousState = savedState;
//  }
}

