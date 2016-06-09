//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4SurfaceList.hh,v 1.6 2006-06-29 18:40:50 gunter Exp $
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
