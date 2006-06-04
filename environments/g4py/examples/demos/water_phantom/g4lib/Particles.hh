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
// $Id: Particles.hh,v 1.2 2006-06-04 21:37:25 kmura Exp $
// ====================================================================
//   Particles.hh
//
//                                         2004 Q
// ====================================================================
#ifndef PARTICLES_H
#define PARTICLES_H

#include "G4VPhysicsConstructor.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class Particles : public G4VPhysicsConstructor {

public:
  Particles();
  ~Particles();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
};

#endif

