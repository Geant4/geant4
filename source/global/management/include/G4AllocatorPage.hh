// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AllocatorPage.hh,v 1.2 1999-11-16 17:40:22 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
