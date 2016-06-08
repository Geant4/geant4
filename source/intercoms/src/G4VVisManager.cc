// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisManager.cc,v 1.2.8.1 1999/12/07 20:49:04 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// Abstract interface for GEANT4 Visualization Manager.
// John Allison 19/Oct/1996.

#include "G4VVisManager.hh"

G4VVisManager::~G4VVisManager () {}

G4VVisManager* G4VVisManager::fpConcreteInstance = 0;
