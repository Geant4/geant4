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
// $Id: G4AssemblyVolume.hh 66872 2013-01-15 01:25:57Z japost $
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
// space.

// Author:      Radovan Chytracek, John Apostolakis, Gabriele Cosmo
// Date:        November 2000
//
// History:
// March 2006, I.Hrivnacova - Extended to support assembly of assemblies
//             of volumes and reflections
// ----------------------------------------------------------------------
#ifndef G4_ASSEMBLYVOLUME_H
#define G4_ASSEMBLYVOLUME_H 

#include <vector>

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
    // The rotation matrix passed as argument can be 0 (identity) or an address
    // even of an object on the upper stack frame. During assembly imprint, a
    // new matrix is created anyway and it is kept track of it so it can be
    // automatically deleted later at the end of the application.
    // This policy is adopted since user has no control on the way the
    // rotations are combined.

  void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                        G4ThreeVector& translation,
                        G4RotationMatrix* rotation);
    //
    // Place the given volume 'pPlacedVolume' inside the assembly.
    //
    // The adopted approach:
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
    // The rotation matrix passed as argument can be 0 (identity) or an address
    // even of an object on the upper stack frame. During assembly imprint, a
    // new matrix is created anyway and it is kept track of it so it can be
    // automatically deleted later at the end of the application.
    // This policy is adopted since user has no control on the way the
    // rotations are combined.

  void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                        G4Transform3D&   transformation);
    //
    // The same as previous, but takes complete 3D transformation in space
    // as its argument.

  void AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                          G4Transform3D&    transformation);
    //
    // The same as previous AddPlacedVolume(), but takes an assembly volume 
    // as its argument.

  void AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                          G4ThreeVector& translation,
                          G4RotationMatrix* rotation);
    //
    // The same as above AddPlacedVolume(), but takes an assembly volume 
    // as its argument with translation and rotation.

  void MakeImprint( G4LogicalVolume* pMotherLV,
                    G4ThreeVector& translationInMother,
                    G4RotationMatrix* pRotationInMother,
                    G4int copyNumBase = 0,
                    G4bool surfCheck = false );
    //
    // Creates instance of an assembly volume inside the given mother volume.

  void MakeImprint( G4LogicalVolume* pMotherLV,
                    G4Transform3D&   transformation,
                    G4int copyNumBase = 0,
                    G4bool surfCheck = false );
    //
    // The same as previous Imprint() method, but takes complete 3D
    // transformation in space as its argument.

  inline std::vector<G4VPhysicalVolume*>::iterator GetVolumesIterator();
  inline unsigned int TotalImprintedVolumes() const;
    //
    // Methods to access the physical volumes imprinted with the assembly.

  unsigned int GetImprintsCount() const;
    //
    // Return the number of made imprints.

  unsigned int GetInstanceCount() const;
    //
    // Return the number of existing instance of G4AssemblyVolume class.

  unsigned int GetAssemblyID()    const;
    //
    // Return instance number of this concrete object.
  
 protected:
     
  void SetInstanceCount( unsigned int value );
  void SetAssemblyID( unsigned int value );
 
  void InstanceCountPlus();
  void InstanceCountMinus();

  void SetImprintsCount( unsigned int value );
  void ImprintsCountPlus();
  void ImprintsCountMinus();
    //
    // Internal counting mechanism, used to compute unique the names of
    // physical volumes created by MakeImprint() methods.

 private:    

  void MakeImprint( G4AssemblyVolume* pAssembly,
                    G4LogicalVolume*  pMotherLV,
                    G4Transform3D&    transformation,
                    G4int copyNumBase = 0,
                    G4bool surfCheck = false );
    //    
    // Function for placement of the given assembly in the given mother
    // (called recursively if the assembly contains an assembly).

 private:

  std::vector<G4AssemblyTriplet> fTriplets;
    //
    // Participating volumes represented as a vector of
    // <logical volume, translation, rotation>.

  std::vector<G4VPhysicalVolume*> fPVStore;
    //
    // We need to keep list of physical volumes created by MakeImprint() method
    // in order to be able to cleanup the objects when not needed anymore.
    // This requires the user to keep assembly objects in memory during the
    // whole job or during the life-time of G4Navigator, logical volume store
    // and physical volume store keep pointers to physical volumes generated by
    // the assembly volume.
    // When an assembly object is about to die it will destroy all its
    // generated physical volumes and rotation matrices as well !

  unsigned int fImprintsCounter;
    //
    // Number of imprints of the given assembly volume.

  static G4ThreadLocal unsigned int fsInstanceCounter;
    //
    // Class instance counter.

  unsigned int fAssemblyID;
    //
    // Assembly object ID derived from instance counter at construction time.

};

#include "G4AssemblyVolume.icc"

#endif // G4_ASSEMBLYVOLUME_H
