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
// $Id: G4Allocator.hh,v 1.11 2002-06-21 16:59:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// A class for fast allocation of objects to the heap through paging
// mechanism. It's meant to be used by associating it to the object to
// be allocated and defining for it new and delete operators via
// MallocSingle() and FreeSingle() methods.
       
//      ---------------- G4Allocator ----------------
//                by Tim Bell, September 1995
// ------------------------------------------------------------

#ifndef G4Allocator_h
#define G4Allocator_h 1

#include <stdlib.h>
#include <stddef.h>

#include "G4AllocatorPage.hh"

template <class Type>
class G4Allocator
{
  public:  // with description

    G4Allocator();
    ~G4Allocator();
      // Constructor & destructor

    inline Type* MallocSingle();
    inline void FreeSingle(Type* anElement);
      // Malloc and Free methods to be used when overloading
      // new and delete operators in the client <Type> object

  private:

    void AddNewPage();
    Type* AddNewElement();	

  private:

    enum { Allocated = 0x47416C, Deleted = 0xB8BE93 };

    G4AllocatorPage<Type> * fPages;
    G4AllocatorUnit<Type> * fFreeList;

    size_t fUnitSize, fPageSize;
};

// ------------------------------------------------------------
// Inline implementation
// ------------------------------------------------------------

// ************************************************************
// G4Allocator constructor
// ************************************************************
//
template <class Type>
G4Allocator<Type>::G4Allocator()
{
  fPages = 0;
  fFreeList = 0;
  fUnitSize = sizeof(G4AllocatorUnit<Type>);
  fPageSize = ( (fUnitSize < 512) ? 1024 : (fUnitSize*10) );
  AddNewPage();
  return;
}

// ************************************************************
// G4Allocator destructor
// ************************************************************
//
template <class Type>
G4Allocator<Type>::~G4Allocator()
{
  G4AllocatorPage<Type> * aPage;
  G4AllocatorPage<Type> * aNextPage;

  aPage = fPages;
  while (aPage != 0)
  {
    aNextPage = aPage->fNext;
    free(aPage->fUnits);
    delete aPage;
    aPage = aNextPage;
  }
  fPages = 0;
  fFreeList = 0;
  return;
}

// ************************************************************
// MallocSingle
// ************************************************************
//
template <class Type>
Type* G4Allocator<Type>::MallocSingle()
{
  Type * anElement;

  if (fFreeList != 0)
  {
    fFreeList->deleted = Allocated;
    anElement = &fFreeList->fElement;
    fFreeList = fFreeList->fNext;
  }
  else
    anElement = AddNewElement();
  return anElement;
}

// ************************************************************
// FreeSingle
// ************************************************************
//
template <class Type>
void G4Allocator<Type>::FreeSingle(Type* anElement)
{
  G4AllocatorUnit<Type> * fUnit;

  // The gcc-3.1 compiler will complain and not correctly handle offsets
  // computed from non-POD types. Pointers to member data should be used
  // instead. This advanced C++ feature seems not to work on earlier
  // versions of the same compiler.
  //
  #if (GNU_GCC==1) && (__GNUC__==3) && (__GNUC_MINOR__>0)
    Type G4AllocatorUnit<Type>::*pOffset = &G4AllocatorUnit<Type>::fElement;
    fUnit = (G4AllocatorUnit<Type> *) ((char *)anElement - size_t(pOffset));
  #else
    fUnit = (G4AllocatorUnit<Type> *)
            ((char *) anElement - offsetof(G4AllocatorUnit<Type>, fElement));
  #endif

  if (fUnit->deleted == Allocated)
  {
    fUnit->deleted = Deleted;
    fUnit->fNext = fFreeList;
    fFreeList = fUnit;
  }
}

// ************************************************************
// AddNewPage
// ************************************************************
//
template <class Type>
void G4Allocator<Type>::AddNewPage()
{
  G4AllocatorPage<Type> * aPage;
  register G4int unit_no;

  aPage = new G4AllocatorPage<Type>;
  aPage->fNext = fPages;
  aPage->fUnits = (G4AllocatorUnit<Type> *)
    malloc(fPageSize);
  fPages = aPage;

  for (unit_no = 0;
       unit_no < G4int(fPageSize/fUnitSize - 1);
       ++unit_no)
  {
    aPage->fUnits[unit_no].fNext = &aPage->fUnits[unit_no + 1];
  }
  aPage->fUnits[unit_no].fNext = fFreeList;
  fFreeList = &aPage->fUnits[0];
}

// ************************************************************
// AddNewElement
// ************************************************************
//
template <class Type>
Type* G4Allocator<Type>::AddNewElement()
{
  Type* anElement;

  AddNewPage();
  fFreeList->deleted = Allocated;
  anElement = &fFreeList->fElement;
  fFreeList=fFreeList->fNext;
  return anElement;
}

#endif
