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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DynamicParticleVector.hh,v 1.5 2001-08-01 17:12:40 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	History: first implementation, alternative to G4FastVector
//               less fast, but it has variable array length and checks boundaries
//	26th September, Chr. Voelcker
// ------------------------------------------------------------

#ifndef G4DynamicParticleVector_h
#define G4DynamicParticleVector_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4DynamicParticle;
#include "g4rw/tpordvec.h"

// #ifdef STL
// //in future use STL vector as container of dynamic particles ...
// typedef Vector<G4DynamicParticle> G4DynamicParticleVector;
// #elseifdef RWT

typedef G4RWTPtrOrderedVector<G4DynamicParticle> G4DynamicParticleVector;

// #endif

#endif
