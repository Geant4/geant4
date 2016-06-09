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
// $Id: G4NavigationLevelRep.cc,v 1.1 2005/06/06 13:41:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//  1 October 1997 J.Apostolakis Initial version. 
//                        
// ----------------------------------------------------------------------

#include "G4NavigationLevelRep.hh"

G4Allocator<G4NavigationLevelRep> aNavigLevelRepAllocator;

// Constructors
//--------------

G4NavigationLevelRep::G4NavigationLevelRep( G4VPhysicalVolume* pPhysVol,
                                      const G4AffineTransform& afTransform,
                                            EVolume            volTp,
                                            G4int              repNo )
   :  sTransform(afTransform),
      sPhysicalVolumePtr(pPhysVol),
      sReplicaNo(repNo),
      sVolumeType(volTp),
      fCountRef(1) 
{
}

G4NavigationLevelRep::G4NavigationLevelRep()
   :  sTransform(),
      sPhysicalVolumePtr(0),
      sReplicaNo(-1),
      sVolumeType(kReplica),
      fCountRef(1) 
{
}

G4NavigationLevelRep::G4NavigationLevelRep( G4VPhysicalVolume* pPhysVol,
                                      const G4AffineTransform& levelAbove,
                                      const G4AffineTransform& relativeCurrent,
                                            EVolume            volTp,
                                            G4int              repNo )
   :  sPhysicalVolumePtr(pPhysVol),
      sReplicaNo(repNo),
      sVolumeType(volTp),
      fCountRef(1) 
{
  sTransform.InverseProduct( levelAbove, relativeCurrent );
}

G4NavigationLevelRep::G4NavigationLevelRep( G4NavigationLevelRep& right )
   :  sTransform(right.sTransform), 
      sPhysicalVolumePtr(right.sPhysicalVolumePtr),
      sReplicaNo(right.sReplicaNo),
      sVolumeType(right.sVolumeType),
      fCountRef(1) 
{
}

// Destructor
//--------------

G4NavigationLevelRep::~G4NavigationLevelRep()
{
#ifdef DEBUG_NAVIG_LEVEL
  if(fCountRef>0)
  {
    G4cerr << "ERROR! - A G4NavigationLevelRep is being deleted that has"
           << G4endl << " positive reference count (fCountRef > 0) !"
           << G4endl;
    G4Exception("G4NavigationLevelRep::~G4NavigationLevelRep()",
                "MemoryCorruption", FatalException,
                "Deletion of data-level object with positive reference count.");
  } 
#endif
}

// Operators
// --------------

G4NavigationLevelRep& 
G4NavigationLevelRep::operator=( const G4NavigationLevelRep &right )
{ 
  if ( &right != this )
  {
    sTransform =  right.sTransform;  
    sPhysicalVolumePtr = right.sPhysicalVolumePtr;
    sVolumeType = right.sVolumeType;
    sReplicaNo =  right.sReplicaNo;
    fCountRef = right.fCountRef;
  }
  return *this;
} 
