// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationLevel.hh,v 1.6 1999-05-17 13:24:06 fbehner Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationLevel
//
// Maintains one level of the geometrical hierarchy. A utility class for use 
//   by G4NavigationHistory.
//
// History:
//
// 30.09.97 J.Apostolakis Initial version. Services derived from
//                        requirements of touchables & G4NavigatorHistory.

#ifndef G4NAVIGATIONLEVEL_HH
#define G4NAVIGATIONLEVEL_HH

// #include <assert.h>
#include "globals.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"

#include "G4NavigationLevelRep.hh"
#include "G4Allocator.hh"

class ostream;

class G4NavigationLevel
{

 public:
   G4NavigationLevel(G4VPhysicalVolume*       newPtrPhysVol,
		     const G4AffineTransform& newT,
		     EVolume                  newVolTp,
		     G4int                    newRepNo= -1);

  // As the previous constructor, but instead of giving Transform, give 
  // the AffineTransform to the level above and the current level's 
  // Transform relative to that.
  //
   G4NavigationLevel(G4VPhysicalVolume*       newPtrPhysVol,
		     const G4AffineTransform& levelAbove,
		     const G4AffineTransform& relativeCurrent,
		     EVolume                  newVolTp,
		     G4int                    newRepNo= -1);

   G4NavigationLevel();
   G4NavigationLevel( const G4NavigationLevel& );

   ~G4NavigationLevel();

   G4NavigationLevel& operator=(const G4NavigationLevel &right);

   G4VPhysicalVolume*       GetPhysicalVolume() const;
   const G4AffineTransform* GetTransformPtr() const ;  // New
   const G4AffineTransform& GetTransform() const ;     // Old
   //    G4AffineTransform& GetTransform();            // Only temporarily
   EVolume                  GetVolumeType() const ;
   G4int                    GetReplicaNo() const ;

   //  To try to resolve the possible problem with returning a reference.
   const G4AffineTransform* GetPtrTransform() const;

   // Override "new" and "delete" to use "G4Allocator".
   inline void *operator new(size_t);
   inline void operator delete(void *aLevel);

   // Pre-allocated 'new' and 'delete' for use with STL 
   //    - do not (directly) use Allocator              F.Behner/J.Apostolakis
   inline void *operator new(size_t, void *);       // void*);
   inline void operator delete(void *ptr, void*);  // Not accepted Sun/HP

 //  Data members: 
 // 
 private:
   G4NavigationLevelRep*  fLevelRep;
};

// extern G4Allocator<G4NavigationLevel> aNavigationLevelAllocator;

#include "G4NavigationLevel.icc"

#endif
