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
// $Id: G4ParticleVector.hh,v 1.6 2001-10-04 20:00:43 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
// HPW decoupling theo models from RW (Mon Mar 16 1998)
// ------------------------------------------------------------

#ifndef G4ParticleVector_h
#define G4ParticleVector_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4DynamicParticle.hh"
#include "g4std/vector"

// #ifdef STL
// for future use STL vector as container 
// typedef Vector<G4DynamicParticle> G4ParticleVector;
// #elseifdef RWT

typedef G4std::vector<G4DynamicParticle *> G4ParticleVector;

#endif
