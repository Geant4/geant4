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

#include "G4VInteractiveSession.hh"

#include "G4InteractorMessenger.hh"

/***************************************************************************/
G4VInteractiveSession::G4VInteractiveSession() { messenger = new G4InteractorMessenger(this); }

/***************************************************************************/
G4VInteractiveSession::~G4VInteractiveSession() { delete messenger; }

/***************************************************************************/
void G4VInteractiveSession::AddMenu(const char*, const char*) {}

/***************************************************************************/
void G4VInteractiveSession::AddButton(const char*, const char*, const char*) {}

/***************************************************************************/
void G4VInteractiveSession::DefaultIcons(bool) {}

/***************************************************************************/
void G4VInteractiveSession::AddIcon(const char*, const char*, const char*, const char*) {}

/***************************************************************************/
void G4VInteractiveSession::SetOutputStyle(const char*, const char*)
{
  G4Exception("G4VInteractiveSession::SetOutputStyle", "uiqt0001", JustWarning,
    "This type of session does not support output styles.");
}

/***************************************************************************/
void G4VInteractiveSession::NativeMenu(bool) {}

/***************************************************************************/
void G4VInteractiveSession::ClearMenu() {}

/***************************************************************************/
void G4VInteractiveSession::UpdateSceneTree(const G4SceneTreeItem&) {}

/***************************************************************************/
void G4VInteractiveSession::AddInteractor(G4String a_name, G4Interactor a_interactor)
{
  interactors[a_name] = a_interactor;
}

/***************************************************************************/
G4Interactor G4VInteractiveSession::GetInteractor(G4String a_name)
{
  G4interactor_map::iterator it;
  if ((it = interactors.find(a_name)) == interactors.end()) return nullptr;
  return (*it).second;
}

/***************************************************************************/
const std::map<G4String, G4VInteractiveSession::OutputStyle>&
G4VInteractiveSession::GetOutputStyles() const
{
  return fOutputStyles;
}

void G4VInteractiveSession::SetStyleUtility(const G4String& destination, const G4String& style)
{
  G4String destinationG4(destination);
  G4String styleG4(style);

  // Lambda expression for changing styles
  const auto& setStyle = [&](const G4String& dest) {
    auto& styleForThisDestination = fOutputStyles[dest];
    if (styleG4 == "fixed") {
      styleForThisDestination.fixed = true;
    }
    else if (styleG4 == "proportional") {
      styleForThisDestination.fixed = false;
    }
    else if (styleG4 == "bold") {
      styleForThisDestination.bold = true;
    }
    else if (styleG4 == "plain") {
      styleForThisDestination.bold = false;
      styleForThisDestination.highlight = false;
    }
    else if (styleG4 == "highlight") {
      styleForThisDestination.highlight = true;
    }
  };

  // Here is where the lambda expression is used
  if (destinationG4 == "all") {
    for (auto& i : fOutputStyles) {
      setStyle(i.first);
    }
  }
  else {
    if (fOutputStyles.find(destinationG4) != fOutputStyles.end()) {
      setStyle(destinationG4);
    }
    else {  // Shouldn't happen, but...
      G4ExceptionDescription ed;
      ed << "Unrecognised output destination \"" << destinationG4 << '"';
      G4Exception("G4VInteractiveSession::SetStyleUtility", "uiqt0002", JustWarning, ed);
    }
  }
}
