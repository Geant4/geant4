// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Allocator.hh,v 1.4 2000-11-20 17:26:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------------- G4Allocator ----------------
//                by Tim Bell, September 1995
// ------------------------------------------------------------
// SG, HPW: Protection vs double deletion of the same element, June 97.

#ifndef G4Allocator_h
#define G4Allocator_h 1

#include <stdlib.h>
#include <stddef.h>

// G4AllocatorPage
#include "G4AllocatorPage.hh"

template <class Type>
class G4Allocator
{
  G4AllocatorPage<Type> *fPages;
  G4AllocatorUnit<Type> *fFreeList;
  
private:
  void AddNewPage();
  Type *AddNewElement();	
     
  enum { Allocated = 0x47416C, Deleted = 0xB8BE93 };

public:
  G4Allocator();
  ~G4Allocator();

  inline Type *MallocSingle()
  {
    Type *anElement;

    if (fFreeList != NULL)
    {
      fFreeList->deleted = Allocated;
      anElement = &fFreeList->fElement;
      fFreeList = fFreeList->fNext;
    }
    else
      anElement = AddNewElement();
    return anElement;
  }

  inline void FreeSingle(Type *anElement)
  {
    G4AllocatorUnit<Type> *fUnit;

    fUnit = (G4AllocatorUnit<Type> *)
      ((char *) anElement -
       offsetof(G4AllocatorUnit<Type>, fElement));
    if (fUnit->deleted == Allocated) {
      fUnit->deleted = Deleted;
      fUnit->fNext = fFreeList;
      fFreeList = fUnit;
    }
/*
    else if  (fUnit->deleted == Deleted) {
      // G4cerr << "G4Allocator : This object is already deleted"  << G4endl;
    } else {
      // G4cerr <<  "G4Allocator: This object is allocated not by G4Allocator"<< G4endl;
    }
*/
  }

};


template <class Type>
G4Allocator<Type>::G4Allocator()
{
  fPages = NULL;
  fFreeList = NULL;
  AddNewPage();
  return;
}

template <class Type>
G4Allocator<Type>::~G4Allocator()
{
  G4AllocatorPage<Type> *aPage;
  G4AllocatorPage<Type> *aNextPage;

  aPage = fPages;
  while (aPage != NULL)
  {
    aNextPage = aPage->fNext;
    free(aPage->fUnits);
    free(aPage);
    aPage = aNextPage;
  }
  fPages = NULL;
  fFreeList = NULL;
  return;
}

static const G4int G4AllocatorPageSize = 1024;

template <class Type>
void G4Allocator<Type>::AddNewPage()
{
  G4AllocatorPage<Type> *aPage;
  register unsigned int unit_no;

  aPage = new G4AllocatorPage<Type>;
  aPage->fNext = fPages;
  aPage->fUnits = (G4AllocatorUnit<Type> *)
    malloc(G4AllocatorPageSize);
  fPages = aPage;

  for (unit_no = 0;
       unit_no < (G4AllocatorPageSize /
		  sizeof(G4AllocatorUnit<Type>)-1);
       ++unit_no)
  {
    aPage->fUnits[unit_no].fNext = &aPage->fUnits[unit_no + 1];
  }
  aPage->fUnits[unit_no].fNext = fFreeList;
  fFreeList = &aPage->fUnits[0];
}

template <class Type>
Type *G4Allocator<Type>::AddNewElement()
{
  Type *anElement;

  AddNewPage();
  fFreeList->deleted = Allocated;
  anElement = &fFreeList->fElement;
  fFreeList=fFreeList->fNext;
  return anElement;
}


#endif
