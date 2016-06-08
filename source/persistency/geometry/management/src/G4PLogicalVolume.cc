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
// $Id: G4PLogicalVolume.cc,v 1.6.2.1 2001/06/28 19:11:28 gunter Exp $
// GEANT4 tag $Name:  $
//
//
//                                      Takashi.Sasaki@kek.jp

#include "G4VSolid.hh"
#include "G4PVSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"
#include "G4PVPhysicalVolume.hh"

#include "G4Material.hh"

#include "HepODBMS/clustering/HepClustering.h"


G4PLogicalVolume::G4PLogicalVolume(
                     const G4LogicalVolume* aLogicalVolume,
                     HepRef(G4PVSolid) persSolid)
 : NoOfDaughters(0)
{
  // Construct persistent logical volume based on the information of
  // the transient logical volume aLogicalVolume.
  // It assumes a persistent persSolid has already been created.
  // Note that daughter volumes will be set by G4PersistencyManager.

  fName = aLogicalVolume->GetName();
    //    *fSensitiveDetector =  LogicalVolume->GetSensitiveDetector();
  fSolid = persSolid; 
    //    *fUserLimits  = LogicalVolume->GetUserLimits();
    //    *fVoxel   = LogicalVolume->GetVoxelHeader();
    //    *fVisAttributes =  LogivalVolume->GetVisAttributes();
    //    *fFastSimulationManager = LogicalVolume->GetFastSimulationManager();
    //      *fIsEnvelops = LogicalVolume->GetIsEnvelope();
  fMaterialName = aLogicalVolume->GetMaterial()->GetName();

}

G4PLogicalVolume::~G4PLogicalVolume()
{
  // ~G4PLogicalVolume() deletes the persistent Logical Volume permanently
  // from the database (as well as from the memory).

  // need to de-reference HepRef's in persistent Physics Volume here.
  // ....

  HepRef(G4PLogicalVolume) aVol = ooThis();

  HepDelete(aVol);
}


//#include "G4MagneticField.hh"
//#include "G4Material.hh"
//#include "G4VSensitiveDetector.hh"
//#include "G4VSolid.hh"
//#include "G4UserLimits.hh"
//#include "G4SmartVoxelHeader.hh"
//#include "G4VisAttributes.hh"
//#include "G4FastSimulationManager.hh"

class G4MagneticField;
class G4Material;
class G4VSensitiveDetector;
class G4VSolid;
class G4UserLimits;
class G4SmartVoxelHeader;
class G4VisAttributes;
class G4FastSimulationManager;

G4LogicalVolume* G4PLogicalVolume::MakeTransientObject(
                                             G4VSolid* theSolid,
                                             G4Material* theMaterial)
{
  // These are null poineters temporaly 
  G4SmartVoxelHeader*    fVoxel = 0;
  const G4VisAttributes* fVisAttributes = 0;

//    G4LogicalVolume(G4VSolid *pSolid, G4Material *pMaterial,
//		    const G4String& name,
//		    G4MagneticField *pField=0,
//		    G4VSensitiveDetector *pSDetector=0,
//		    G4UserLimits *pULimits=0);

  G4String name;
  G4LogicalVolume* aLogicalVolume = new G4LogicalVolume(
                       theSolid, theMaterial, name=fName, 0, 0, 0);

  // if ( aLogicalVolume )
  //   aLogicalVume->SetVisAttributes(*fVisAttributes);

  return aLogicalVolume;
}
//////////////////////////

G4int G4PLogicalVolume::GetNoDaughters() const
{  
//	return fDaughters.entries();
  return NoOfDaughters;
}

HepRef(G4PVPhysicalVolume) G4PLogicalVolume::GetDaughter(const G4int i) const
{  
  return fDaughters[i];
}

void G4PLogicalVolume::AddDaughter(HepRef(G4PVPhysicalVolume) p)
{
  fDaughters.insert_element(p);
  NoOfDaughters++; 
}

G4bool G4PLogicalVolume::IsDaughter(const HepRef(G4PVPhysicalVolume) p) const
{
  for ( G4int i=0; i<NoOfDaughters; i++) 
    if ( fDaughters[i] == p )
      return true;

  return false;
}

void G4PLogicalVolume::RemoveDaughter(const HepRef(G4PVPhysicalVolume) p)
{
  for ( G4int i=0; i<NoOfDaughters; i++) 
    if ( fDaughters[i] == p )
    {
      fDaughters.replace_element_at(0,i);
      NoOfDaughters--; 
    }
}

HepRef(G4PVSolid) G4PLogicalVolume::GetSolid()
{
  return fSolid;
}

void G4PLogicalVolume::SetSolid(HepRef(G4PVSolid) pSolid)
{
  fSolid = pSolid;
}
