// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceList.hh,v 1.1 1999-01-07 16:07:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4SurfaceList_h
#define __G4SurfaceList_h 1

#include "G4Surface.hh"

class G4SurfaceList
{
 public:
  
  G4SurfaceList();
  ~G4SurfaceList();
  
  int number_of_elements;
  
  G4Surface* first;
  G4Surface* next;
  G4Surface* last;
  G4Surface* temp;
  G4Surface* index;
  
  void MoveToFirst(G4Surface *srf);
  void AddSurface(G4Surface *srf);
  
  G4Surface* GetSurface();
  G4Surface* GetSurface(int number);
  G4Surface* GetLastSurface();
  
  void RemoveSurface(G4Surface* srf);
  void RemovePointer();
  
  void MoveToFirst();
  void Step();
  
  void EmptyList();
  void G4SortList();
  
  void QuickG4Sort(G4Surface**, int, int);

};

#endif








