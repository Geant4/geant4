// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSVolumeHandle.hh,v 1.1 2001-03-16 16:47:57 radoone Exp $
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
// For more details see G4ReferenceCountedHandle.hh file.
//
// Author:      Radovan Chytracek
// Version:     1.0
// Date:        March 2001
// ----------------------------------------------------------------------
#ifndef _G4GRSVOLUMEHANDLE_H_
#define _G4GRSVOLUMEHANDLE_H_ 1

#include "G4GRSVolume.hh"
#include "G4ReferenceCountedHandle.hh"

typedef G4ReferenceCountedHandle<G4GRSVolume> G4GRSVolumeHandle;

#endif // _G4GRSVOLUMEHANDLE_H_
