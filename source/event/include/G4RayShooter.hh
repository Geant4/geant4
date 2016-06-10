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
//
// $Id: G4RayShooter.hh 66892 2013-01-17 10:57:59Z gunter $
//

// class description:
//
//  the object of this class shoots a ray (actually geantino) primary particle
// associated with a G4Event object. This slass must be used exclusively by
// G4RayTracer.
//

#ifndef G4RayShooter_h
#define G4RayShooter_h 1


#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;

class G4RayShooter
{
  public: 
     G4RayShooter();
  public:
     virtual ~G4RayShooter();

  public: // with description
     void Shoot(G4Event* evt,G4ThreeVector vtx,G4ThreeVector direc);
     // This method generates a primary vertex and a primary particle at the
     // given vertex point and with the given direction. This method is invoked
     // by G4RayTracer and G4MaterialScanner.

  private:  
     void SetInitialValues();

     G4ParticleDefinition* particle_definition;
     G4ParticleMomentum    particle_momentum_direction;
     G4double	           particle_energy;
     G4ThreeVector         particle_position;
     G4double              particle_time;
     G4ThreeVector         particle_polarization;
};

#endif







