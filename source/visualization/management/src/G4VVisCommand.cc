// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisCommand.cc,v 1.2 1999-01-09 16:31:17 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

G4VisManager* G4VVisCommand::fpVisManager = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneNotifyHandlers = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneRemove = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneSelect = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneHandlerAttach = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneHandlerRemove = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandSceneHandlerSelect = 0;

G4UIcommand*        G4VVisCommand::fpCommandViewerCreate = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandViewerRemove = 0;

G4UIcmdWithAString* G4VVisCommand::fpCommandViewerSelect = 0;
