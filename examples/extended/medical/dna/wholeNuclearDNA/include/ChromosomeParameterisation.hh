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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file ChromosomeParameterisation.hh
/// \brief Definition of the ChromosomeParameterisation class


#ifndef ChromosomeParameterisation_H
#define ChromosomeParameterisation_H 1

#include <G4VPVParameterisation.hh>
#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4Torus.hh>
#include <G4Trap.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>
#include <G4Ellipsoid.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include <vector>

class G4VPhysicalVolume;
class G4Box;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;
class G4Ellipsoid;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChromosomeParameterisation : public G4VPVParameterisation
{
public:
  ChromosomeParameterisation(const char* filename);

  virtual ~ChromosomeParameterisation();
  inline int GetNumRosettes()
  {
    return fPositions.size();
  }

  virtual void
  ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

  virtual void
  ComputeDimensions(G4Tubs& rosette,
      const G4int copyNo,
      const G4VPhysicalVolume* physVol) const;

private:
  // Dummy declarations to get rid of warnings ...
  virtual void
  ComputeDimensions(G4Box&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Trd&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Trap&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Sphere&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Orb&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Torus&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Para&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Hype&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Polycone&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Polyhedra&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Ellipsoid&, const G4int, const G4VPhysicalVolume*) const
  {
  }

  std::vector<G4ThreeVector*> fPositions;
  std::vector<G4RotationMatrix*> fRotations;
};

#endif
