// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: rw2stl_algo.h,v 1.4 1999/11/25 10:14:45 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RW2STL_LessPtr
//
//  STL wrapper utility class G4RW2STL_LessPtr.
//---------------------------------------------------------------

#ifndef rw2stl_algo_h
#define rw2stl_algo_h

#include "g4rw/defs.h"

template <class T>
class G4RW2STL_LessPtr
{
public:

  G4bool operator()(const T* a, const T* b) const
    { 
      if(a==0) return false;
      if(b==0) return true;
      return *a<*b;
    }

};

#endif
