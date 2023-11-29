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
// G4StateManager class implementation
//
// Authors: G.Cosmo, M.Asai - November 1996
// --------------------------------------------------------------------

#include "G4StateManager.hh"
#include "G4ios.hh"

// Initialization of the static pointer of the single class instance
//
G4ThreadLocal G4StateManager* G4StateManager::theStateManager = nullptr;
G4int G4StateManager::verboseLevel                            = 0;

// --------------------------------------------------------------------
G4StateManager::G4StateManager()
{
#ifdef G4MULTITHREADED
  G4iosInitialization();
#endif
}

G4StateManager::~G4StateManager()
{
  G4VStateDependent* state = nullptr;

  while(!theDependentsList.empty())
  {
    state = theDependentsList.back();
    theDependentsList.pop_back();
    for(auto i = theDependentsList.cbegin(); i != theDependentsList.cend();)
    {
      if(*i == state)
      {
        i = theDependentsList.erase(i);
      }
      else
      {
        ++i;
      }
    }
    delete state;
  }
  theStateManager = nullptr;
#ifdef G4MULTITHREADED_DEACTIVATE
  G4iosFinalization();
#endif
}

// --------------------------------------------------------------------
G4StateManager* G4StateManager::GetStateManager()
{
  if(theStateManager == nullptr)
  {
    theStateManager = new G4StateManager;
  }
  return theStateManager;
}

// --------------------------------------------------------------------
G4bool G4StateManager::RegisterDependent(G4VStateDependent* aDependent,
                                         G4bool bottom)
{
  G4bool ack = true;
  if(!bottom)
  {
    theDependentsList.push_back(aDependent);
  }
  else
  {
    if(theBottomDependent != nullptr)
    {
      theDependentsList.push_back(theBottomDependent);
    }
    theBottomDependent = aDependent;
  }
  return ack;
}

// --------------------------------------------------------------------
G4bool G4StateManager::DeregisterDependent(G4VStateDependent* aDependent)
{
  G4VStateDependent* tmp = nullptr;
  for(auto i = theDependentsList.cbegin(); i != theDependentsList.cend();)
  {
    if(**i == *aDependent)
    {
      tmp = *i;
      i   = theDependentsList.erase(i);
    }
    else
    {
      ++i;
    }
  }
  return (tmp != nullptr);
}

// --------------------------------------------------------------------
const G4ApplicationState& G4StateManager::GetCurrentState() const
{
  return theCurrentState;
}

// --------------------------------------------------------------------
const G4ApplicationState& G4StateManager::GetPreviousState() const
{
  return thePreviousState;
}

// --------------------------------------------------------------------
G4bool G4StateManager::SetNewState(const G4ApplicationState& requestedState)
{
  return SetNewState(requestedState, nullptr);
}

// --------------------------------------------------------------------
G4bool G4StateManager::SetNewState(const G4ApplicationState& requestedState,
                                   const char* msg)
{
  if(requestedState == G4State_Abort && suppressAbortion > 0)
  {
    if(suppressAbortion == 2)
    {
      return false;
    }
    if(theCurrentState == G4State_EventProc)
    {
      return false;
    }
  }
  msgptr                        = msg;
  std::size_t i                 = 0;
  G4bool ack                    = true;
  G4ApplicationState savedState = thePreviousState;
  thePreviousState              = theCurrentState;

  while((ack) && (i < theDependentsList.size()))
  {
    ack = theDependentsList[i]->Notify(requestedState);
    ++i;
  }
  if(theBottomDependent != nullptr)
  {
    ack = theBottomDependent->Notify(requestedState);
  }

  if(!ack)
  {
    thePreviousState = savedState;
  }
  else
  {
    theCurrentState = requestedState;
    if(verboseLevel > 0)
    {
      G4cout << "#### G4StateManager::SetNewState from "
             << GetStateString(thePreviousState) << " to "
             << GetStateString(requestedState) << G4endl;
    }
  }
  msgptr = nullptr;
  return ack;
}

// --------------------------------------------------------------------
G4VStateDependent* G4StateManager::RemoveDependent(
  const G4VStateDependent* aDependent)
{
  G4VStateDependent* tmp = nullptr;
  for(auto i = theDependentsList.cbegin(); i != theDependentsList.cend();)
  {
    if(**i == *aDependent)
    {
      tmp = *i;
      i   = theDependentsList.erase(i);
    }
    else
    {
      ++i;
    }
  }
  return tmp;
}

// --------------------------------------------------------------------
G4String G4StateManager::GetStateString(const G4ApplicationState& aState) const
{
  G4String stateName;
  switch(aState)
  {
    case G4State_PreInit:
      stateName = "PreInit";
      break;
    case G4State_Init:
      stateName = "Init";
      break;
    case G4State_Idle:
      stateName = "Idle";
      break;
    case G4State_GeomClosed:
      stateName = "GeomClosed";
      break;
    case G4State_EventProc:
      stateName = "EventProc";
      break;
    case G4State_Quit:
      stateName = "Quit";
      break;
    case G4State_Abort:
      stateName = "Abort";
      break;
    default:
      stateName = "Unknown";
      break;
  }
  return stateName;
}

// --------------------------------------------------------------------
void G4StateManager::SetVerboseLevel(G4int val) { verboseLevel = val; }
