// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolume.cc,v 1.3 1999-12-15 14:49:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4LogicalVolume Implementation
//
// History:
// 12.02.99 S.Giani: Deafult initialization of voxelization quality.
// 10.20.97 - P. MoraDeFreitas : "Fast" replaces "Parameterisation" in
//            class/method names. (release B.00 for parameterisation).
// 04.08.97 P.MoraDeFreitas/J.A. Added methods for ParameterisedSimulation 
// 19.08.96 P.Kent Modified for G4VSensitive Detector
// 11.07.95 P.Kent Initial version

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"

// Constructor - set member data and add to logical Store, zero voxel ptr
//               Initialises daughter vector to 0 length
G4LogicalVolume::G4LogicalVolume( G4VSolid *pSolid,  G4Material *pMaterial,
				  const G4String& name,
				  G4FieldManager *pFieldMgr,
				  G4VSensitiveDetector *pSDetector,
				  G4UserLimits *pULimits) : 
   fDaughters(0), fVoxel(0), fSmartless(2.), fVisAttributes (0) , fFastSimulationManager (0),
   fIsEnvelope(FALSE), fFieldManager(pFieldMgr)
{
    SetSolid(pSolid);
    SetMaterial(pMaterial);
    SetName(name);
    SetSensitiveDetector(pSDetector);
    SetUserLimits(pULimits);    
    // Add to solid Store
    G4LogicalVolumeStore::Register(this);
}

// Destructor - remove from solid Store
// NOTE: Not virtual
G4LogicalVolume::~G4LogicalVolume()
{
    G4LogicalVolumeStore::DeRegister(this);
}

//  As this method is recursive, inlining it is harder (and asking for it 
//  is pointless if not counterproductive: it will increase code size)
//
void G4LogicalVolume::SetFastSimulationManager (
	G4FastSimulationManager* pNewFastSimul,
	G4bool IsEnvelope) 
{
  if(!fIsEnvelope || IsEnvelope) {
     fIsEnvelope=IsEnvelope;
     fFastSimulationManager = pNewFastSimul;

     G4int NoDaughters=GetNoDaughters();
     while((NoDaughters--)>0){
        G4LogicalVolume *DaughterLogVol; 
        DaughterLogVol= GetDaughter(NoDaughters)->GetLogicalVolume();
        if(DaughterLogVol->GetFastSimulationManager() != pNewFastSimul) {
           DaughterLogVol->SetFastSimulationManager(pNewFastSimul,FALSE);
        }
     }
  }
}

void G4LogicalVolume::
ClearEnvelopeForFastSimulation( G4LogicalVolume* motherLogVol) 
{
  if( fIsEnvelope ) {
     G4FastSimulationManager* NewFastSimulationVal=NULL;

     // This is no longer an envelope !
     fIsEnvelope=FALSE;

     if( motherLogVol == NULL ) {
	motherLogVol = this->FindMotherLogicalVolumeForEnvelope();

     } else {
        //  Check that motherLogVol is this' mother.
        //   If not, raise exception and set it to NULL.
        if ( FALSE ){
	   G4Exception(
            "G4LogicalVolume::ClearEnvelope Gave wrong mother LogicalVolume" );
	   motherLogVol = NULL;
	}
     }
     // Reset the ParameterisedSimulation values of self and all daughters
     //    (after ensuring the mother was given correctly or was found)
     if( motherLogVol != NULL ) {
	  NewFastSimulationVal = motherLogVol->
	    GetFastSimulationManager();
	  this->SetFastSimulationManager (NewFastSimulationVal,
					    FALSE);
     }
  }else{
     G4Exception("Called G4LogicalVolume::ClearEnvelope for a non-envelope Logical Volume" );
  }
}

void G4LogicalVolume::SetFieldManager (
	G4FieldManager* pNewFieldMgr,
        G4bool          forceAllDaughters) 
{
  fFieldManager = pNewFieldMgr;

  G4int NoDaughters=GetNoDaughters();
  while((NoDaughters--)>0){
    G4LogicalVolume *DaughterLogVol; 
    DaughterLogVol= GetDaughter(NoDaughters)->GetLogicalVolume();
    if(forceAllDaughters || (DaughterLogVol->GetFieldManager() != 0)) {
      DaughterLogVol->SetFieldManager(pNewFieldMgr, forceAllDaughters);
    }
  }

}

//  The following method returns a meaningful result IF and only IF
//  the current logical volume has exactly one physical volume that
//  uses it.

G4LogicalVolume* 
G4LogicalVolume::FindMotherLogicalVolumeForEnvelope()
{
  G4LogicalVolume* motherLogVol= 0;
  G4LogicalVolumeStore *Store = G4LogicalVolumeStore::GetInstance();

  // Look for the current volume's mother volume.

  for (G4int LV=0;LV < Store->entries(); LV++){
     G4LogicalVolume *aLogVol;
    
     aLogVol= Store->at(LV);
     if( (aLogVol!=this) &&         // Don't look for it inside itself...
	(aLogVol->GetFastSimulationManager()!=NULL)){
      
        for (G4int daughter=0; daughter < aLogVol->GetNoDaughters();
	     daughter++){
	   if( aLogVol->GetDaughter(daughter)->GetLogicalVolume()==this) { 
	      // "- Oh Dear, aLogVol is my mother !!!"
	      motherLogVol = aLogVol;
	      break;
	   }
	}
     }
  }
  return motherLogVol;
}
