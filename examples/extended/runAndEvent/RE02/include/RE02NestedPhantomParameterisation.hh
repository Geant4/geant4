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
// $Id: RE02NestedPhantomParameterisation.hh,v 1.4 2010-11-09 19:46:30 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VNestedParameterisation
//
// Class description:
//
// Base class for parameterisations that use information from the parent
// volume to compute the material of a copy/instance of this volume. 
// This is in addition to using the current replication number.
// 
// Notes:
//  - Such a volume can be nested inside a placement volume or a parameterised  
//    volume.
//  - The user can modify the solid type, size or transformation using only
//    the replication number of this parameterised volume.
//    He/she is NOT allowed to change these attributes using information of
//    parent volumes - otherwise incorrect results will occur.
//  Also note that the usual restrictions apply: 
//   - the mother volume, in which these copies are placed, must always be
//     of the same dimensions

// History:
// 24.02.05 - J.Apostolakis - First created version.
// --------------------------------------------------------------------
#ifndef RE02NESTEDPARAMETERISATION_HH
#define RE02NESTEDPARAMETERISATION_HH

#include "G4Types.hh"
#include "G4VNestedParameterisation.hh" 
#include "G4ThreeVector.hh"
#include <vector>

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

// CSG Entities which may be parameterised/replicated
//
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

class RE02NestedPhantomParameterisation: public G4VNestedParameterisation
{
  public:  // with description

    RE02NestedPhantomParameterisation(const G4ThreeVector& voxelSize,
				      G4int nz,
				      std::vector<G4Material*>& mat);
    virtual ~RE02NestedPhantomParameterisation(); 

    // Methods required in derived classes
    // -----------------------------------
    virtual G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
					const G4int repNo, 
                                        const G4VTouchable *parentTouch=0
                                        );
  // Required method, as it is the reason for this class.
  //   Must cope with parentTouch=0 for navigator's SetupHierarchy

    virtual G4int       GetNumberOfMaterials() const;
    virtual G4Material* GetMaterial(G4int idx) const;
      // Needed to define materials for instances of Nested Parameterisation 
      //   Current convention: each call should return the materials 
      //   of all instances with the same mother/ancestor volume.

    virtual void ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume *currentPV) const;

    // Methods optional in derived classes
    // -----------------------------------

    // Additional standard Parameterisation methods, 
    //   which can be optionally defined, in case solid is used.

    virtual void ComputeDimensions(G4Box &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const;

private:  // Dummy declarations to get rid of warnings ...
  void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) 
    const {}
  void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*)
    const {}
  void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) 
    const {}
//  G4Material* ComputeMaterial(const G4int repNo,
//                              G4VPhysicalVolume* currentVol,
//                              const G4VTouchable* parentTouch)
//  { return ComputeMaterial( currentVol, repNo, parentTouch ); }
  using G4VNestedParameterisation::ComputeMaterial;

private:
  G4double fdX,fdY,fdZ;
  G4int fNz;
  //
  std::vector<G4double>  fpZ;
  std::vector<G4Material*> fmat;
};

#endif
