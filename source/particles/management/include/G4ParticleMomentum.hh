// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleMomentum.hh,v 1.3 1999-10-28 23:24:12 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ------------------------------------------------------------
//
// Class Description
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  G4ParticleMomentum is "momentum direction" not "momentum vector"
//  The name is miss-leading so you should not use G4ParticleMomentum
//  and you are recommended to use G4ThreeVector instead
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//

#ifndef G4ParticleMomentum_h
#define G4ParticleMomentum_h 1

#include "G4ThreeVector.hh"

typedef G4ThreeVector G4ParticleMomentum;
//  Is a three vector containing the momentum direction.
//  The momentum components are in global reference system.

#endif






