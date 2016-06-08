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
// $Id: tvvector.h,v 1.8.4.1 2001/06/28 19:09:59 gunter Exp $
// GEANT4 tag $Name:  $
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
