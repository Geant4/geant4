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
// $Id: tpsrtvec.h,v 1.8.8.1 2001/06/28 19:09:58 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWTPtrSortedVector
//
//  Class description:
//
//  STL wrapper class for Pointer Sorted arrays. It implements
//  Rogue Wave RWTPtrSortedVector signature but intrinsically
//  using STL vector.

//---------------------------------------------------------------

#ifndef __tpsrtvec
#define __tpsrtvec

#include "g4rw/rw2stl_algo.h"
#include "g4std/vector"
#include "g4std/algorithm"
#include "g4rw/defs.h"


template <class T>
class G4RWTPtrSortedVector : public G4std::vector<T*>
{
 
  typedef G4std::vector<T*> std_pvector;
  typedef typename std_pvector::iterator iterator;
  typedef typename std_pvector::const_iterator const_iterator;
 
public:

  G4RWTPtrSortedVector(size_t n=G4RWDEFAULT_CAPACITY);
  G4RWTPtrSortedVector(const G4RWTPtrSortedVector<T>&);
  virtual ~G4RWTPtrSortedVector();

  inline G4RWTPtrSortedVector<T>& operator=(const G4RWTPtrSortedVector<T>&);

  inline T* const& operator [] (size_t) const;
  inline T* const& operator () (size_t) const;

  inline T*& at (size_t);
  inline T*  at (size_t) const;

  void clear();
  void clearAndDestroy ();

  T* find (const T*) const;
  inline size_t entries () const;
  inline void insert (T*);
  size_t index (const T* a);
  inline G4bool isEmpty () const;

  inline T* const & first () const; 
  inline T* const & last () const;

  size_t occurrencesOf(const T*) const;

  T* remove (const T*);
  size_t removeAll (const T*);

private:

  size_t rwsize; 
  G4RW2STL_LessPtr<T> sorter;  
	
};

// Include actual implementations of templated class methods.

#include "g4rw/tpsrtvec.icc"

#endif
