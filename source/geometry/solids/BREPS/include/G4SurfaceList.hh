// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceList.hh,v 1.3 2000-08-28 08:57:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SurfaceList
//
// Class description:
// 
// Class defining a list of surfaces.

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef __G4SurfaceList_h
#define __G4SurfaceList_h 1

#include "G4Surface.hh"

class G4SurfaceList
{
 public:  // with description
  
  G4SurfaceList();
  ~G4SurfaceList();
    // Constructor & destructor.

  void MoveToFirst(G4Surface *srf);
  void AddSurface(G4Surface *srf);
  
  G4Surface* GetSurface();
  G4Surface* GetSurface(G4int number);
  G4Surface* GetLastSurface();
  
  void RemoveSurface(G4Surface* srf);
  void RemovePointer();
  
  void MoveToFirst();
  void Step();
  
  void EmptyList();
  void G4SortList();
  
  void QuickG4Sort(G4Surface**, G4int, G4int);

 public:  // without description

  G4int number_of_elements;
  
  G4Surface* first;
  G4Surface* next;
  G4Surface* last;
  G4Surface* temp;
  G4Surface* index;
  
};

#endif
