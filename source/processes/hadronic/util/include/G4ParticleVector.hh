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
// $Id: G4ParticleVector.hh,v 1.8 2002-12-12 19:18:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ParticleVector_h
#define G4ParticleVector_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4DynamicParticle.hh"
#include "g4std/vector"

typedef G4std::vector<G4DynamicParticle *> G4ParticleVector;
struct DeleteDynamicParticle {void operator()(G4DynamicParticle*aP) {delete aP;}};

#endif
