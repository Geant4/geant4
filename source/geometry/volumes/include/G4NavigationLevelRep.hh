// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationLevelRep.hh,v 1.4 2000-04-25 16:15:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationLevelRep
//
// Class description:
//
// A data representation class, used to hold the data for a single level
// of the Navigation history tree.
//
// This is the body of a handle/body pair of classes, that implement
// reference counting for NavigationLevels.
// The corresponding handle class is G4NavigationLevel

// History:
//
// - 1 October 1997, J.Apostolakis: initial version. 

#ifndef G4NAVIGATIONLEVELREP_HH
#define G4NAVIGATIONLEVELREP_HH

#include "globals.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Allocator.hh"

class G4NavigationLevelRep
{

 public:  // with description

   G4NavigationLevelRep(G4VPhysicalVolume*       newPtrPhysVol,
		        const G4AffineTransform& newT,
		        EVolume                  newVolTp,
		        G4int                    newRepNo= -1);

   G4NavigationLevelRep(G4VPhysicalVolume*       newPtrPhysVol,
		        const G4AffineTransform& levelAbove,
		        const G4AffineTransform& relativeCurrent,
		        EVolume                  newVolTp,
		        G4int                    newRepNo= -1);
     // As the previous constructor, but instead of giving Transform, give 
     // the AffineTransform to the level above and the current level's 
     // Transform relative to that.

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

   void   AddAReference(); 
   G4bool RemoveAReference(); 
     // Take care of the reference counts.

   inline void *operator new(size_t);
     // Override "new"    to use "G4Allocator".
   inline void operator delete(void *aTrack);
     // Override "delete" to use "G4Allocator".

 public:  // without description

   // const G4AffineTransform* GetPtrTransform() const;
   // G4AffineTransform&  accessTransform();

 private:

   G4AffineTransform  sTransform;
     // Compounded global->local transformation (takes a point in the 
     // global reference system to the system of the volume at this level)

   G4VPhysicalVolume* sPhysicalVolumePtr;
     // Physical volume ptrs, for this level's volume

   G4int              sReplicaNo;
   EVolume            sVolumeType;
     // Volume `type' 

   G4int              fCountRef; 

};

#include "G4NavigationLevelRep.icc"

// extern G4Allocator<G4NavigationLevelRep> aNavigLevelRepAllocator;

#endif
