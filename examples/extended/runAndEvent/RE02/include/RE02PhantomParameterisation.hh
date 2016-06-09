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
// $Id: RE02PhantomParameterisation.hh,v 1.2 2006/06/29 17:44:58 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#ifndef RE02PhantomParameterisation_h
#define RE02PhantomParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4VisAttributes;
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

//
class RE02PhantomParameterisation : public G4VPVParameterisation
{
public:

  RE02PhantomParameterisation(const G4ThreeVector& motherSize,
			     const G4int nx, const G4int ny, const G4int nz); 

  virtual ~RE02PhantomParameterisation();

  void ComputeTransformation(const G4int copyNo, 
			     G4VPhysicalVolume* physVol)const;

  private:  // Dummy declarations to get rid of warnings ...
  void ComputeDimensions(G4Box&, const G4int,const G4VPhysicalVolume* ) const{};
  void ComputeDimensions(G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Para&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

private:
  G4ThreeVector fDxyzMother;
  G4ThreeVector fDxyz;
  G4int         fNx, fNy, fNz;
  std::vector<G4ThreeVector> fPositions;

};
#endif

