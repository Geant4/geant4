// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationLevel.hh,v 1.12 2000-11-01 16:51:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationLevel
//
// Class description:
//
// Maintains one level of the geometrical hierarchy.
// A utility class for use by G4NavigationHistory.

// History:
//
// 30.09.97 J.Apostolakis Initial version. Services derived from
//                        requirements of touchables & G4NavigatorHistory.

#ifndef G4NAVIGATIONLEVEL_HH
#define G4NAVIGATIONLEVEL_HH

#include "globals.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"

#include "G4NavigationLevelRep.hh"
#include "G4Allocator.hh"

class G4NavigationLevel
{

 public:  // with description

   G4NavigationLevel(G4VPhysicalVolume*       newPtrPhysVol,
		     const G4AffineTransform& newT,
		     EVolume                  newVolTp,
		     G4int                    newRepNo= -1);

   G4NavigationLevel(G4VPhysicalVolume*       newPtrPhysVol,
		     const G4AffineTransform& levelAbove,
		     const G4AffineTransform& relativeCurrent,
		     EVolume                  newVolTp,
		     G4int                    newRepNo= -1);
     // As the previous constructor, but instead of giving Transform, give 
     // the AffineTransform to the level above and the current level's 
     // Transform relative to that.

   G4NavigationLevel();
   G4NavigationLevel( const G4NavigationLevel& );

   ~G4NavigationLevel();

   G4NavigationLevel& operator=(const G4NavigationLevel &right);

   inline G4VPhysicalVolume*       GetPhysicalVolume() const;
   inline const G4AffineTransform* GetTransformPtr() const ;  // New
   inline const G4AffineTransform& GetTransform() const ;     // Old

   inline EVolume                  GetVolumeType() const ;
   inline G4int                    GetReplicaNo() const ;

 public:  // without description

   // inline    G4AffineTransform& GetTransform();

   inline const G4AffineTransform* GetPtrTransform() const;
     // To try to resolve the possible problem with returning a reference.

   inline void *operator new(size_t);
   inline void operator delete(void *aLevel);
     // Override "new" and "delete" to use "G4Allocator".

   inline void *operator new(size_t, void *);
#ifndef G4NOT_ISO_DELETES
   inline void operator delete(void *ptr, void*);  // Not accepted Sun/HP
#endif
     // Pre-allocated 'new' and 'delete' for use with STL 
     // Do not (directly) use Allocator   -   F.Behner/J.Apostolakis

 private:

   G4NavigationLevelRep*  fLevelRep;
};

// extern G4Allocator<G4NavigationLevel> aNavigationLevelAllocator;

#include "G4NavigationLevel.icc"

#endif
