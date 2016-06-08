// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPersistentSubDbMan.cc,v 1.2 1999/11/26 20:37:00 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// class G4VPersistentSubDbMan 
//
// History:
// 99.11.25 Y.Morita  Initial version

#include "G4VPersistentSubDbMan.hh"

G4VPersistentSubDbMan::G4VPersistentSubDbMan()
 : f_DB(NULL), f_container(NULL), f_currentContainer(NULL)
{;}

G4VPersistentSubDbMan::~G4VPersistentSubDbMan()
{;}

