// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TouchableHandle.hh,v 1.1 2001-03-13 16:06:23 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Class G4TouchableHandle
//
// Class description:
//
// A class to provide reference counting mechanism for any kind of touchable
// objects.
// The basic rule for the use of this class is that it is passed always by
// reference and is not constructed by calling new.
//
// Before the use of this smart pointer object one can test its validity
// using the operator !().
//  if( !smartPtrObj ) { ... } // OK!
//  else               { ... } // Problem! We must initialize it first!
//
// The code which tries to delete this object won't compile, because it is
// not a pointer it is an object.
//
// Author:      Radovan Chytracek
// Version:     1.0
// Date:        February 2001
// ----------------------------------------------------------------------
#ifndef _G4TOUCHABLEHANDLE_H_
#define _G4TOUCHABLEHANDLE_H_ 1

#include "G4VTouchable.hh>
#include "G4ReferenceCountedHandle.hh"

typedef G4ReferenceCountedHandle<G4VTouchable> G4TouchableHandle;

#endif // _G4TOUCHABLEHANDLE_H_
