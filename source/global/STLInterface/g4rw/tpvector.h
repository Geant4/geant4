// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tpvector.h,v 1.9 2000/02/11 14:09:24 gcosmo Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWTPtrVector
//
//  Class description:
//
//  STL wrapper class for Pointer arrays. It implements
//  Rogue Wave RWTPtrVector signature but intrinsically
//  using STL vector.

//---------------------------------------------------------------

#ifndef __tpvector
#define __tpvector

#include "g4std/vector"
#include "g4rw/defs.h"

template<class T>
class G4RWTPtrVector : public G4std::vector<T*>
{
 
  typedef G4std::vector<T*> std_pvector;
  typedef typename std_pvector::iterator iterator;
  typedef typename std_pvector::const_iterator const_iterator;
 
public:

  G4RWTPtrVector ();
  G4RWTPtrVector (size_t);
  G4RWTPtrVector (size_t, T* const &);
  G4RWTPtrVector (const G4RWTPtrVector<T>&);
  virtual ~G4RWTPtrVector();

  G4RWTPtrVector<T>& operator = (T*);
  
  inline T* operator () (size_t) const;
  inline T*& operator () (size_t);

  inline T* operator [] (size_t) const;
  inline T*& operator [] (size_t);
  inline T* const * data() const;
  inline size_t length();
  inline void reshape(size_t);
  void resize(size_t);

private:

  size_t rwsize;

};

// Include actual implementations of templated class methods.

#include "g4rw/tpvector.icc"

#endif
