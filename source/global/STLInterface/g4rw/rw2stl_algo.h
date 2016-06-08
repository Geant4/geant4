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
// $Id: rw2stl_algo.h,v 1.4.8.1 2001/06/28 19:09:58 gunter Exp $
// GEANT4 tag $Name:  $
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
