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
// $Id: G4ReferenceCountedHandle.cc,v 1.1 2001-11-06 17:09:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4ReferenceCountedHandle.cc
//
//  Author:
//    18 Nov 2001, Gabriele Cosmo - CERN
// --------------------------------------------------------------

#include "G4ReferenceCountedHandle.hh"

template <class X>
G4Allocator<G4ReferenceCountedHandle<X>::CountedObject >
  G4ReferenceCountedHandle<X>::CountedObject::aCountedObjectAllocator;

template <class X>
G4Allocator<G4ReferenceCountedHandle<X> >
  G4ReferenceCountedHandle<X>::aRCHAllocator;

