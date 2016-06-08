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
// $Id: tpordvec.h,v 1.7.8.1 2001/06/28 19:09:58 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWTPtrOrderedVector
//
//  Class description:
//
//  STL wrapper class for Pointer Ordered arrays. It implements
//  Rogue Wave RWTPtrOrderedVector signature but intrinsically
//  using STL vector.

//---------------------------------------------------------------

#ifndef __tpordvec
#define __tpordvec

#include "g4std/vector"
#include "g4std/set"
#include "g4std/functional"
#include "g4rw/defs.h"

template<class T>
class G4RWTPtrOrderedVector : public G4std::vector<T*>
{
  typedef G4std::vector<T*> std_pvector;
  typedef typename std_pvector::iterator iterator;
  typedef typename std_pvector::const_iterator const_iterator;

public:

  G4RWTPtrOrderedVector(size_t capacity=G4RWDEFAULT_CAPACITY);
  G4RWTPtrOrderedVector(const G4RWTPtrOrderedVector<T>&);
  virtual ~G4RWTPtrOrderedVector();

  inline G4RWTPtrOrderedVector<T>& operator=(const G4RWTPtrOrderedVector<T>&);

  inline T*& operator()(size_t);
  inline T* const& operator()(size_t) const;
  inline T*& operator[](size_t);
  inline T* const& operator[](size_t) const;

  void clearAndDestroy();
  void clear();

  size_t index (const T*) const;
  void resize (size_t);
  inline void insert (T*);
  void insertAt (size_t, T*); 
  size_t occurrencesOf (const T*) const;	

  T* remove (const T*);	
  size_t removeAll (const T*);	
  T* removeAt (size_t);
  T* removeFirst ();
  T* removeLast ();

  T* find (const T*) const;															
  inline size_t entries() const;
  inline void append (T*);
  inline size_t length() const;
  inline T*& at (size_t);
  inline T*const& at (size_t) const;
  G4bool contains (const T*) const;
  inline G4bool isEmpty() const;
  
  inline T* first() const;
  inline T* last() const;
  inline void prepend(T* const);

private:

  size_t rwsize;

};

// Include actual implementations of templated class methods.

#include "g4rw/tpordvec.icc"

#endif
