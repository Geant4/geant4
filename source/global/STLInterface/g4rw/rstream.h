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
// $Id: rstream.h,v 1.4.8.1 2001/06/28 19:09:58 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
//---------------------------------------------------------------
//  GEANT 4 utility file
//
//  Definition of utility function rwEatwhite.
//---------------------------------------------------------------

#ifndef __rstream
#define __rstream

#include "g4rw/defs.h"
#include "G4Types.hh"
#include "g4std/iostream"

inline G4std::istream& rwEatwhite(G4std::istream& stream)
{
  return stream >> G4std::ws;
}

#endif
