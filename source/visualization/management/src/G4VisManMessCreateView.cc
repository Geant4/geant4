// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessCreateView.cc,v 1.2 1999-01-09 16:31:29 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Create View sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"
#include <string.h>

void G4VisManMessenger::AddCommandCreateView () {

  G4UIcommand* command;
  G4UIparameter* param;

  ///////////////////  /vis~/create_view/new_graphics_system  /////
  //cr_vw \hline
  //cr_vw /vis~/create\_view/ new\_graphics\_system & choice &
  //cr_vw Creates a new scene and a new view of a new graphics system; all
  //cr_vw become current. \\%
  command = new G4UIcommand    ("/vis~/create_view/new_graphics_system", this);
  command -> SetGuidance
    (
     "Creates a new scene and a new view of a new graphics system; all "
     "become current."
     );
  param   =  new G4UIparameter ("Graphics system selector", 's', true);
  param   -> SetDefaultValue  ("none");
  const G4GraphicsSystemList& gslist =
    fpVMan -> GetAvailableGraphicsSystems ();
  G4String candidates;
  for (int igslist = 0; igslist < gslist.entries (); igslist++) {
    G4String nickname = gslist (igslist) -> GetNickname ();
    if (nickname . length () > 0) {
      candidates += nickname;
    }
    else {
      candidates += gslist (igslist) -> GetName ();
    }
    candidates += " ";
  }
  param   -> SetParameterCandidates(candidates);
  command -> SetParameter (param);
  fCommandList.append (command);

}

void G4VisManMessenger::DoCommandCreateView (const G4String& commandPath,
					     G4String& newValues) {

  ///////////////////  /vis~/create_view/new_graphics_system  /////
  if (commandPath == "/vis~/create_view/new_graphics_system") {
    G4String selector;
    const char* aString = newValues;
    istrstream is ((char*) aString) ; is >> selector;
    const G4GraphicsSystemList& gsl = fpVMan -> GetAvailableGraphicsSystems ();
    int nSystems = gsl.entries ();
    if (nSystems > 0) {
      int iGS;  // Selector index.
      for (iGS = 0; iGS < nSystems; iGS++) {
	if (selector.compareTo (gsl [iGS] -> GetName (),
				RWCString::ignoreCase) == 0 ||
	    selector.compareTo (gsl [iGS] -> GetNickname (),
				RWCString::ignoreCase) == 0) {
	  break;  // Match found.
	}
      }
      if (iGS < 0 || iGS >= nSystems) {
	// Invalid command line argument or non.  Print available systems.
	G4bool help;
	if (selector.compareTo ("help", RWCString::ignoreCase) == 0 ||
	    selector.compareTo ("info", RWCString::ignoreCase) == 0)
	  help = true;
	else help = false;
	int maxNameLength = strlen ("class name");
	int nameLength;
	int maxNicknameLength = strlen ("nickname");
	int nicknameLength;
	int igs;
	for (igs = 0; igs < nSystems; igs++) {
	  nameLength = strlen (gsl [igs] -> GetName ());
	  if (maxNameLength < nameLength) maxNameLength = nameLength;
	  nicknameLength = strlen (gsl [igs] -> GetNickname ());
	  if (maxNicknameLength < nicknameLength)
	    maxNicknameLength = nicknameLength;
	}
	G4String className, nickname, description;
	className = "class name";
	int i;
	for (i = className.length (); i < maxNameLength; i++)
	  className = className + " ";
	G4cout << "Available graphics systems:\n  " << className
	     << "nickname";
	className = "Help";
	for (i = className.length (); i < maxNameLength; i++)
	  className = className + " ";
	if (help) G4cout << '\n';
	G4cout << "\n " << className << "  info        for more information";
	for (igs = 0; igs < nSystems; igs++) {
	  className = gsl [igs] -> GetName ();
	  for (i = className.length (); i < maxNameLength; i++)
	    className = className + " ";
	  if (help) G4cout << '\n';
	  G4cout << "\n " << className;
	  nickname = gsl [igs] -> GetNickname ();
	  for (i = nickname.length (); i < maxNameLength; i++)
	    nickname = nickname + " ";
	  G4cout << "  " << nickname;
	  if (help) {
	    description = gsl [igs] -> GetDescription ();
	    if (description != "") {
	      size_t index = 0;
	      while ((index = description.index ('\n', ++index)) != RW_NPOS) {
		description.replace (index, 1, "\n   ");
	      }
	      G4cout << "\n   " << description;
	    }
	  }
	}
	if (nSystems == 1) iGS = 0;  // Only one system - create view anyway.
	else {
	  G4cout <<
	    "\nChoose by specifying class name or nickname, if any, or help.";
	}
	G4cout << endl;
      }
      if (iGS >=0 && iGS < nSystems) {
	// Valid index.  Create view.
	G4VGraphicsSystem* pSystem = gsl [iGS];;
	fpVMan -> SetCurrentGraphicsSystemAndCreateView (pSystem);
	if (fpVMan -> GetVerboseLevel () > 0) {
	  G4cout << "Graphics system set to " << pSystem -> GetName () << endl;
	  if (fpVMan -> GetVerboseLevel () > 1) {
	    fpVMan -> PrintCurrentSystem ();
	  }
	}
      }
    }
  }
}
