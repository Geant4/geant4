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
// $Id: G4NavigationLevelRep.hh 85846 2014-11-05 15:45:28Z gcosmo $
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
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONLEVELREP_HH
#define G4NAVIGATIONLEVELREP_HH

#include "G4Types.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Allocator.hh"

#include "geomwdefs.hh"

class G4NavigationLevelRep
{

 public:  // with description

   inline G4NavigationLevelRep( G4VPhysicalVolume*  newPtrPhysVol,
                          const G4AffineTransform&  newT,
                                EVolume             newVolTp,
                                G4int               newRepNo= -1 );

   inline G4NavigationLevelRep( G4VPhysicalVolume*  newPtrPhysVol,
                          const G4AffineTransform&  levelAbove,
                          const G4AffineTransform&  relativeCurrent,
                                EVolume             newVolTp,
                                G4int               newRepNo= -1 );
     // As the previous constructor, but instead of giving Transform, give 
     // the AffineTransform to the level above and the current level's 
     // Transform relative to that.

   inline G4NavigationLevelRep();
   inline G4NavigationLevelRep( G4NavigationLevelRep& );

   inline ~G4NavigationLevelRep();

   inline G4NavigationLevelRep& operator=(const G4NavigationLevelRep &right);

   inline G4VPhysicalVolume* GetPhysicalVolume();

   inline const G4AffineTransform* GetTransformPtr() const ;  // New
   inline const G4AffineTransform& GetTransform() const ;     // Old

   inline EVolume            GetVolumeType() const ;
   inline G4int              GetReplicaNo() const ;

   inline void   AddAReference(); 
   inline G4bool RemoveAReference(); 
     // Take care of the reference counts.

   inline void *operator new(size_t);
     // Override "new"    to use "G4Allocator".
   inline void operator delete(void *aTrack);
     // Override "delete" to use "G4Allocator".

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

#endif
