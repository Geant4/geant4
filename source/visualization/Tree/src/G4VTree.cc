// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTree.cc,v 1.1 2001-04-10 15:09:35 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#include "G4VTree.hh"
#include "G4VTreeSceneHandler.hh"
#include "G4VTreeViewer.hh"

G4VTree::G4VTree (const G4String& name,
		  const G4String& nickname,
		  const G4String& description,
		  Functionality f):
  G4VGraphicsSystem (name,
		     nickname,
		     description,
		     f) {}

G4VTree::~G4VTree () {}
