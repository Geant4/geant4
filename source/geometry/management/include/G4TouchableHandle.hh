//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TouchableHandle.hh,v 1.5 2001-11-06 17:06:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Class G4TouchableHandle
//
// Class description:
//
// A type providing reference counting mechanism for any kind of touchable
// objects.
// The basic rule for the use of this type is that this handle must always
// be exchanged by reference never dinamically allocated (i.e. never
// instantiated using 'new').
// For more details see G4ReferenceCountedHandle.

// Author:      Radovan Chytracek, CERN  (Radovan.Chytracek@cern.ch)
// Version:     1.0
// Date:        February 2001
// ----------------------------------------------------------------------
#ifndef _G4TOUCHABLEHANDLE_H_
#define _G4TOUCHABLEHANDLE_H_ 1

#include "G4VTouchable.hh"
#include "G4ReferenceCountedHandle.hh"

typedef G4ReferenceCountedHandle<G4VTouchable> G4TouchableHandle;

#endif // _G4TOUCHABLEHANDLE_H_
