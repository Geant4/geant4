// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisCommand.cc,v 1.6 1999-12-16 17:19:31 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

G4VVisCommand::~G4VVisCommand () {}

G4VisManager* G4VVisCommand::fpVisManager = 0;

G4std::vector <G4UIcommand*> G4VVisCommand::sceneNameCommands;

G4std::vector <G4UIcommand*> G4VVisCommand::sceneHandlerNameCommands;

G4std::vector <G4UIcommand*> G4VVisCommand::viewerNameCommands;
