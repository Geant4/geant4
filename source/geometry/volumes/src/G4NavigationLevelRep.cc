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
// $Id: G4NavigationLevelRep.cc,v 1.2 2006-06-29 18:58:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
