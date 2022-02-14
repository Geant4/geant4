//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRParallelWorldPhysics.hh
//   Utility class that adds GRParallelWorldBiasingProcess to the
//   process managers of particles
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRParallelWorldPhysics_h
#define GRParallelWorldPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class GRParallelWorldPhysics : public G4VPhysicsConstructor
{
public:

  GRParallelWorldPhysics(const G4String& name = "ParallelWP", G4bool layerdMass = false);
  virtual ~GRParallelWorldPhysics();

public:

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

   // hide assignment operator
  GRParallelWorldPhysics & operator=(const GRParallelWorldPhysics &right);
  GRParallelWorldPhysics(const GRParallelWorldPhysics&);

  G4bool fLayeredMass;
};

#endif
