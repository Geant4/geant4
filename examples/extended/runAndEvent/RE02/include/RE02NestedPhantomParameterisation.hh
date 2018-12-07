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
/// \file runAndEvent/RE02/include/RE02NestedPhantomParameterisation.hh
/// \brief Definition of the RE02NestedPhantomParameterisation class
//
//
//
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
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

//
/// A nested parameterisation class for a phantom
///
///  (Description)
///     This parameterisation handles material and transfomation of voxles.
///
/// - G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
///                               const G4int repNo,
///                               const G4VTouchable *parentTouch=0)
///     returns material.
///       if ix%2==0 && iy%2==0 && iz%2==0 then fMat[0]
///       else fMat[1]
///
/// - G4int GetNumberOfMaterials() const
///     returns the number of material defined in fMat
///
/// - G4Material* GetMaterial(G4int idx) const
///     returns the i-th material of fMat
///
/// - void ComputeTransformation(const G4int no,
///                              G4VPhysicalVolume *currentPV) const
///     returns a transformation with the physical volume of the 2nd argument
///     according to copyNo
///     Its position is defined as G4ThreeVector(0.,0.,fpZ[copyNo]).
///
/// - void ComputeDimensions(G4Box &, const G4int, 
///                          const G4VPhysicalVolume *) const
///     returns dimensions of this parameterized volume with the physical 
///     volume of the 3rd argument.
///
//
class RE02NestedPhantomParameterisation: public G4VNestedParameterisation
{
  public:  // with description

    RE02NestedPhantomParameterisation(const G4ThreeVector& voxelSize,
                                      G4int nz,
                                      std::vector<G4Material*>& mat);
   ~RE02NestedPhantomParameterisation(); 

    // Methods required in derived classes
    // -----------------------------------
    G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
                                const G4int repNo, 
                                const G4VTouchable *parentTouch=0
                                        );
  // Required method, as it is the reason for this class.
  //   Must cope with parentTouch=0 for navigator's SetupHierarchy

    G4int       GetNumberOfMaterials() const;
    G4Material* GetMaterial(G4int idx) const;
      // Needed to define materials for instances of Nested Parameterisation 
      //   Current convention: each call should return the materials 
      //   of all instances with the same mother/ancestor volume.

    void ComputeTransformation(const G4int no,
                               G4VPhysicalVolume *currentPV) const;

    // Methods optional in derived classes
    // -----------------------------------

    // Additional standard Parameterisation methods, 
    //   which can be optionally defined, in case solid is used.

    void ComputeDimensions(G4Box &,
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
  void ComputeDimensions (G4Ellipsoid&,const G4int,const G4VPhysicalVolume*) 
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
  std::vector<G4Material*> fMat;
};

#endif
