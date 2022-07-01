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
// G4RayShooter
//
// Class description:
//
// The object of this class shoots a ray (actually geantino) primary particle
// associated with a G4Event object. This slass must be used exclusively by
// G4RayTracer.

// Author: Makoto Asai, 2000
// --------------------------------------------------------------------
#ifndef G4RayShooter_hh
#define G4RayShooter_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;

class G4RayShooter
{
  public:

    G4RayShooter() = default;
    virtual ~G4RayShooter() = default;

    void Shoot(G4Event* evt, G4ThreeVector vtx, G4ThreeVector direc);
      // Generates a primary vertex and a primary particle at the given
      // vertex point and with the given direction. This method is invoked
      // by G4RayTracer and G4MaterialScanner.

  private:
    G4ParticleDefinition* particle_definition = nullptr;
    G4ParticleMomentum    particle_momentum_direction;
    G4double              particle_energy = 1.0*CLHEP::GeV;
    G4ThreeVector         particle_position;
    G4double              particle_time = 0.0;
    G4ThreeVector         particle_polarization;
};

#endif
