// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tvvector.h,v 1.7 1999/11/26 17:17:58 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWTValVector
//
//  Class description:
//
//  STL wrapper class for Value arrays. It implements
//  Rogue Wave RWTValVector signature but intrinsically
//  using STL vector.

//---------------------------------------------------------------

#ifndef __tvvector
#define __tvvector

#include "g4std/vector"
#include "g4rw/defs.h"

template<class T>
class G4RWTValVector : public G4std::vector<T>
{
 
  typedef G4std::vector<T> std_vector;
  typedef typename std_vector::iterator iterator;
  typedef typename std_vector::const_iterator const_iterator;
 
public:

  G4RWTValVector ();
  G4RWTValVector (size_t);
  G4RWTValVector (size_t, const T&);
  G4RWTValVector (const G4RWTValVector<T>&);
  virtual ~G4RWTValVector();

  inline G4RWTValVector<T>& operator= (const G4RWTValVector<T>&);

  inline T& operator()(size_t);
  inline const T& operator()(size_t) const;
  inline T& operator[](size_t);
  inline const T& operator[](size_t) const;

  inline size_t length() const;

  inline const T* data() const;
  void resize(size_t);
  inline void reshape(size_t);
  inline const T& ref(size_t) const;

private:

  size_t rwsize;

};

// Include actual implementations of templated class methods.

#include "g4rw/tvvector.icc"

#endif
