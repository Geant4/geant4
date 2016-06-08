// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tvordvec.h,v 1.8 2000/03/10 15:53:24 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWTValOrderedVector
//
//  Class description:
//
//  STL wrapper class for Value Ordered arrays. It implements
//  Rogue Wave RWTValOrderedVector signature but intrinsically
//  using STL vector.

//---------------------------------------------------------------

#ifndef __tvordvec
#define __tvordvec

#include "g4std/vector"
#include "g4std/algorithm"

#include "g4rw/defs.h"
#include "g4rw/tvvector.h"

template<class T>
class G4RWTValOrderedVector : public G4std::vector<T> 
{
 
  typedef G4std::vector<T> std_vector;
  typedef typename std_vector::iterator iterator;
  typedef typename std_vector::const_iterator const_iterator;
 
public:

  G4RWTValOrderedVector(size_t capacity=G4RWDEFAULT_CAPACITY);
  G4RWTValOrderedVector(const G4RWTValOrderedVector<T>&);
  virtual ~G4RWTValOrderedVector();

  inline G4RWTValOrderedVector<T>& operator=(const G4RWTValOrderedVector<T>&);

  inline T& operator()(size_t);
  inline const T& operator()(size_t) const;
  inline T& operator[](size_t);
  inline const T& operator[](size_t) const;

  inline T& at (size_t);
  inline T  at (size_t) const;

  void resize (size_t);
  void clear( void );
  G4bool contains(const T) const;
  inline size_t length() const;
  size_t index (const T&);

  inline G4bool insert (const T&); 
  void insertAt (size_t, const T&); 

  size_t occurrencesOf (const T&) const;	

  G4bool remove (const T&);	
  size_t removeAll (const T&);	
  T removeAt (size_t);
	
  size_t entries () const;
  inline void append (const T&);

  inline T first() const;
  inline T last() const;

private:

  size_t rwsize;
  
};

// Include actual implementations of templated class methods.

#include "g4rw/tvordvec.icc"

#endif

