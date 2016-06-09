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
// $Id: RE02PhantomParameterisation.hh,v 1.1 2005/11/24 01:44:18 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

