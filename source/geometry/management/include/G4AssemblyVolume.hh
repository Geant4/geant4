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
// $Id: G4AssemblyVolume.hh,v 1.8 2002-09-10 16:59:44 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4AssemblyVolume
//
// Class description:
//
// G4AssemblyVolume is a helper class to make the build process of geometry
// easier. It allows to combine several volumes together in an arbitrary way
// in 3D space and then work with the result as with a single logical volume
// for placement.
// The resulting objects are independent copies of each of the assembled
// logical volumes. The placements are not, however, bound one to each other
// when placement is done. They are seen as independent physical volumes in
// space. This is the main difference when compared to the boolean solids or
// parametrised physical volumes.

// Author:      Radovan Chytracek, John Apostolakis, Gabriele Cosmo
// Version:     1.0
// Date:        November 2000
// ----------------------------------------------------------------------

#ifndef G4_ASSEMBLYVOLUME_H
#define G4_ASSEMBLYVOLUME_H 

#include "g4std/vector"
#include "G4Transform3D.hh"
#include "G4AssemblyTriplet.hh"

class G4VPhysicalVolume;

class G4AssemblyVolume
{
 public:  // with description

  G4AssemblyVolume();    
  G4AssemblyVolume( G4LogicalVolume* volume,
                    G4ThreeVector& translation,
                    G4RotationMatrix* rotation);
  ~G4AssemblyVolume();
    //
    // Constructors & destructor.
    // At destruction all the generated physical volumes and associated
    // rotation matrices of the imprints will be destroyed.
    //
    // The rotation matrix passed in can be 0 = identity or an address even of an object
    // on the upper stack frame. During assembly imprint, it creates anyway a new matrix
    // and keeps track of it so it can delete it later at destruction time.
    // This new policy has been adopted since user has no control on the way the rotations
    // are combined it's safer doing it this way.
    //
    // WARNING! This interface will likely change in the next major release of Geant4 from
    //          a pointer to a reference due to the reason above
    //

  void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                        G4ThreeVector& translation,
                        G4RotationMatrix* rotation);
    //
    // Place the given volume 'pPlacedVolume' inside the assembly.
    //
    // The taken approach:
    //
    // - Place it w.r.t. the assembly coordinate system.
    //   This step is applied to each of the participating volumes.
    //
    // The other possible approaches:
    //
    // - Place w.r.t. the firstly added volume.
    //   When placed the first, the virtual coordinate system becomes
    //   the coordinate system of the first one.
    //   Every next volume being added into the assembly will be placed
    //   w.r.t to the first one.
    //
    // - Place w.r.t the last placed volume.
    //   When placed the first, the virtual coordinate system becomes
    //   the coordinate system of the first one.
    //   Every next volume being added into the assembly will be placed
    //   w.r.t to the previous one.
    //
    // The rotation matrix passed in can be 0 = identity or an address even of an object
    // on the upper stack frame. During assembly imprint, it creates anyway a new matrix
    // and keeps track of it so it can delete it later at destruction time.
    // This new policy has been adopted since user has no control on the way the rotations
    // are combined it's safer doing it this way.
    //
    // WARNING! This interface will likely change in the next major release of Geant4 from
    //          a pointer to a reference due to the reason above
    //

  void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                        G4Transform3D&   transformation);
    //
    // The same as previous but takes complete 3D transformation in space
    // as its argument.

  void MakeImprint( G4LogicalVolume* pMotherLV,
                    G4ThreeVector& translationInMother,
                    G4RotationMatrix* pRotationInMother,
                    G4int copyNumBase = 0 );
    //
    // Creates instance of an assembly volume inside the given mother volume.

  void MakeImprint( G4LogicalVolume* pMotherLV,
                    G4Transform3D&   transformation,
                    G4int copyNumBase = 0 );
    //
    // The same as previous but takes complete 3D transformation in space
    // as its argument.

 private:    

  G4std::vector<G4AssemblyTriplet> fTriplets;
    //
    // Participating volumes represented as a vector of
    // <logical volume, translation, rotation>.

  G4std::vector<G4VPhysicalVolume*> fPVStore;
    //
    // We need to keep list of physical volumes created by MakeImprint() method
    // in order to be able to cleanup the objects when not needed anymore.
    // This requires the user to keep assembly objects in memory during the
    // whole job or during the life-time of G4Navigator, log. vol store
    // and phys. vol. store may keep pointers to physical volumes generated by
    // the assembly volume.
    // When an assembly object is about to die it will destroy all its
    // generated physical volumes and rotation matrices as well !
    // This may affect validity of detector contruction !

 public:

  unsigned int GetImprintsCount() const;
  
 protected:
    
  void         SetImprintsCount( unsigned int value );
  void         ImprintsCountPlus();
  void         ImprintsCountMinus();
    //
    // Internal counting mechanism, used to compute unique the names of
    // phys. volumes created by MakeImprint(...) method(s).

 private:

  unsigned int fImprintsCounter;
    //
    // Number of imprints of the given assembly volume.

 public:
    
  unsigned int GetInstanceCount() const;
    // Return the number of existing instance of G4AssemblyVolume class
  unsigned int GetAssemblyID()    const;
    // Return instance number of this concrete object
  
 protected:
     
  void         SetInstanceCount( unsigned int value );
  void         SetAssemblyID( unsigned int value );
 
  void         InstanceCountPlus();
  void         InstanceCountMinus();

 private:

  static unsigned int fsInstanceCounter;
  // Class instance counter
         unsigned int fAssemblyID;
  // Assembly object ID derived from instance counter at construction time

};

#include "G4AssemblyVolume.icc"

#endif // G4_ASSEMBLYVOLUME_H

