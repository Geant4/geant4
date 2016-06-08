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
// $Id: G4AllocatorPage.hh,v 1.2.8.1 2001/06/28 19:10:01 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------------- G4AllocatorPage ----------------
//                by Tim Bell, September 1995
// ------------------------------------------------------------

#ifndef G4AllocatorPage_h
#define G4AllocatorPage_h 1

#include "G4AllocatorUnit.hh"

template <class Type>
class G4AllocatorPage
{
public:
  G4AllocatorPage<Type> *fNext;
  G4AllocatorUnit<Type> *fUnits;
};

#endif
