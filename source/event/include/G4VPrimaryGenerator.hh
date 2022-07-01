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
// G4VPrimaryGenerator
//
// Class description:
//
// This is an abstract base class of all of primary generators.
// This class has only one pure virtual method GeneratePrimaryVertex()
// which takes a G4Event object and generates a primay vertex and
// primary particles associate to the vertex

// Author: Makoto Asai (SLAC)
// --------------------------------------------------------------------
#ifndef G4VPrimaryGenerator_hh
#define G4VPrimaryGenerator_hh 1

#include "G4ThreeVector.hh"

class G4Event;

class G4VPrimaryGenerator
{
  public:

    G4VPrimaryGenerator() = default;
    virtual ~G4VPrimaryGenerator() = default;
      // Constructor and destructor

    static G4bool CheckVertexInsideWorld(const G4ThreeVector& pos);
      // Static service method for checking a point is included
      // in the (current) world

    virtual void GeneratePrimaryVertex(G4Event* evt) = 0;
      // Pure virtual method which a concrete class derived from this
      // base class must implement

    inline G4ThreeVector GetParticlePosition() { return particle_position; }
    inline G4double GetParticleTime() { return particle_time; }

    inline void SetParticlePosition(G4ThreeVector aPosition)
      { particle_position = aPosition; }
    inline void SetParticleTime(G4double aTime)
      { particle_time = aTime; }

  protected:

    G4ThreeVector particle_position;
    G4double particle_time = 0.0;
};

#endif
