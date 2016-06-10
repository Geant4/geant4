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
// $Id: G4NavigationLevel.hh 86527 2014-11-13 15:06:24Z gcosmo $
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
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONLEVEL_HH
#define G4NAVIGATIONLEVEL_HH

#include "G4Types.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"

#include "G4NavigationLevelRep.hh"
#include "G4Allocator.hh"

#include "geomwdefs.hh"

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
     // Do not (directly) use Allocator

 private:

   G4NavigationLevelRep*  fLevelRep;
};

#include "G4NavigationLevel.icc"

#endif
