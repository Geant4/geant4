// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationLevelRep.hh,v 1.3 1999-12-15 16:40:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationLevelRep
//
//    A data representation class, used to hold the data for
//    a single level of the Navigation history tree.
//
//    This is the body of a handle/body pair of classes,
//     that implement reference counting for NavigationLevels.
//    The corresponding handle class is G4NavigationLevel
//
// History:
//
//  1 October 1997,   J.Apostolakis:    Initial version. 

#ifndef G4NAVIGATIONLEVELREP_HH
#define G4NAVIGATIONLEVELREP_HH

#include "globals.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Allocator.hh"

class G4NavigationLevelRep
{

 public:
   G4NavigationLevelRep(G4VPhysicalVolume*       newPtrPhysVol,
		        const G4AffineTransform& newT,
		        EVolume                  newVolTp,
		        G4int                    newRepNo= -1);

  // As the previous constructor, but instead of giving Transform, give 
  // the AffineTransform to the level above and the current level's 
  // Transform relative to that.
  //
   G4NavigationLevelRep(G4VPhysicalVolume*       newPtrPhysVol,
		        const G4AffineTransform& levelAbove,
		        const G4AffineTransform& relativeCurrent,
		        EVolume                  newVolTp,
		        G4int                    newRepNo= -1);

   G4NavigationLevelRep();
   G4NavigationLevelRep( G4NavigationLevelRep& );

   ~G4NavigationLevelRep();

   G4NavigationLevelRep& operator=(const G4NavigationLevelRep &right);

   G4VPhysicalVolume* GetPhysicalVolume() const;

   const G4AffineTransform* GetTransformPtr() const ;  // New
   const G4AffineTransform& GetTransform() const ;     // Old
   //    G4AffineTransform& GetTransform();            // Only temporarily
   EVolume                  GetVolumeType() const ;
   G4int              GetReplicaNo() const ;

   // const G4AffineTransform* GetPtrTransform() const;
   // G4AffineTransform&  accessTransform();

   //  Take care of the reference counts.
   //
   void   AddAReference(); 
   G4bool RemoveAReference(); 

   inline void *operator new(size_t);
      // Override "new"    to use "G4Allocator".
   inline void operator delete(void *aTrack);
      // Override "delete" to use "G4Allocator".

 private:
   // Compounded global->local transformation (takes a point in the 
   // global reference system to the system of the volume at this level)
   G4AffineTransform  sTransform;
   // Physical volume ptrs, for this level's volume
   G4VPhysicalVolume* sPhysicalVolumePtr;
   // Volume `type' 
   G4int              sReplicaNo;
   EVolume            sVolumeType;

   G4int              fCountRef; 

};

#include "G4NavigationLevelRep.icc"

// extern G4Allocator<G4NavigationLevelRep> aNavigLevelRepAllocator;

#endif
