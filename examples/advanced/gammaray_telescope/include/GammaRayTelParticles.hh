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
// $Id: RemsimParticles.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Remsim example
// ----------------------------------------------------------------------------
// Code developed by:
//
// S. Guatelli(b)
// 
// ----------------------------------------------------------------------------

#ifndef GAMMARAYTELPARTICLES_HH
#define GAMMARAYTELPARTICLES_HH 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class GammaRayTelParticles : public G4VPhysicsConstructor {
public: 

  GammaRayTelParticles(const G4String& name = "particles");
  
  virtual ~GammaRayTelParticles();
  
  virtual void ConstructParticle();
  
  // This method is dummy
  virtual void ConstructProcess() {};  
};
#endif








