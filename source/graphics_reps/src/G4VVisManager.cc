// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisManager.cc,v 1.1 1999-01-07 16:09:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Visualization Manager.
// John Allison 19/Oct/1996.

#include "G4VVisManager.hh"

G4VVisManager* G4VVisManager::fpConcreteInstance = 0;
