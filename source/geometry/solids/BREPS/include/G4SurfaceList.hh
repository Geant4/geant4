// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceList.hh,v 1.4 2000-11-08 14:22:04 gcosmo Exp $
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
  const G4Surface* GetSurface(G4int number);
  const G4Surface* GetLastSurface() const;
  
  void RemoveSurface(G4Surface* srf);
  void RemovePointer();
  
  void MoveToFirst();
  void Step();
  
  void EmptyList();
  void G4SortList();
  
  void QuickG4Sort(G4Surface**, G4int, G4int);

  const G4Surface* GetFirst() const { return first; }
  const G4Surface* GetNext()  const { return next;  }
  G4int GetSize() const { return number_of_elements; }

 private:
 
  G4SurfaceList(const G4SurfaceList&);
  G4SurfaceList& operator=(const G4SurfaceList&);
    // Private copy constructor and assignment operator.

 private:  // without description

  G4int number_of_elements;
  
  G4Surface* first;
  G4Surface* next;
  G4Surface* last;
  G4Surface* temp;
  G4Surface* index;
  
};

#endif
