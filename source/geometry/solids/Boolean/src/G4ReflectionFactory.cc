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
// $Id: G4ReflectionFactory.cc,v 1.2 2001-11-08 15:47:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Ivana Hrivnacova, 16.10.2001  (Ivana.Hrivnacova@cern.ch)

// 
// Class G4ReflectionFactory Implementation
//
// Decomposition of a general transformation
// that can include reflection in a "reflection-free" transformation:
// 
// x(inM') = TG*x(inM)         TG - general transformation
//         = T*(R*x(inM))      T  - "reflection-free" transformation
//         = T* x(inReflM)   
//
// Daughters transformation:
// When a volume V containing daughter D with transformation TD
// is placed in mother M with a general tranformation TGV,
// the TGV is decomposed,
// new reflected volume ReflV containing a new daughter ReflD
// with reflected transformation ReflTD is created:
// 
// x(inV) = TD * x(inD);
// x(inM) = TGV * x(inV) 
//        = TV * R * x(inV) 
//	  = TV * R * TD * x(inD)
//	  = TV * R*TD*R-1 * R*x(inD)
//	  = TV * ReflTD * x(inReflD)


#include "G4ReflectionFactory.hh"
#include "G4ReflectedSolid.hh"  
#include "G4LogicalVolume.hh"  
#include "G4PVPlacement.hh"  
#include "G4PVReplica.hh"  

G4ReflectionFactory* G4ReflectionFactory::fInstance = 0;
const G4String  G4ReflectionFactory::fNameExtension = "_refl";
const G4Scale3D G4ReflectionFactory::fScale = G4ScaleZ3D(-1.0);

//_____________________________________________________________________________
G4ReflectionFactory* G4ReflectionFactory::Instance() 
{
// Static singleton access method.
// ---

  if (!fInstance) new G4ReflectionFactory();
  
  return fInstance;
}  

//_____________________________________________________________________________
G4ReflectionFactory::G4ReflectionFactory()
  : fVerboseLevel(0)     
{
// Protected singleton constructor.
// ---

  fInstance = this;
}

//_____________________________________________________________________________
G4ReflectionFactory::~G4ReflectionFactory() {
}

//
// public methods
//

//_____________________________________________________________________________
G4PhysicalVolumesPair  G4ReflectionFactory::Place(
                                          const G4Transform3D& transform3D,
                                          const G4String&      name,
	   	                          G4LogicalVolume* LV,
					  G4LogicalVolume* motherLV,
					  G4bool  isMany, 
					  G4int   copyNo)
{
// Evaluates the passed transformation; if it contains reflection
// it performs its decomposition, creates new reflected solid and
// logical volume (or retrieves them from a map if the reflected
// objects were already created), transforms the daughters (if present)
// and place it in the given mother.
// The result is a pair of physical volumes;
// the second physical volume is a placement in a reflected mother
// - or 0 if mother LV was not reflected.
// ---

  if (fVerboseLevel>0) {
    G4cout << "Place " << name << " lv " << LV << " " 
           << LV->GetName() << G4endl;
  }  

  // decompose transformation
  G4Scale3D     scale;
  G4Rotate3D    rotation;
  G4Translate3D translation;

  transform3D.getDecomposition(scale, rotation, translation);
  G4Transform3D pureTransform3D = translation * rotation;
  
  //PrintTransform(transform3D);
  //PrintTransform(pureTransform3D);

  // check that scale correspond to fScale
  CheckScale(scale);
  
  //
  // reflection IS NOT present in transform3D 
  //

  if (! IsReflection(scale)) {
    if (fVerboseLevel>0)
      G4cout << "scale positive" << G4endl;

    G4VPhysicalVolume* pv1
      =  new G4PVPlacement(pureTransform3D, LV, name, motherLV, isMany, copyNo);
 
    G4VPhysicalVolume* pv2 = 0;
    if (G4LogicalVolume* reflMotherLV = GetReflectedLV(motherLV)) {

      // if mother was reflected
      // reflect this LV and place it in reflected mother
      
      pv2 = new G4PVPlacement(
                  fScale * (pureTransform3D * fScale.inverse()),
	          ReflectLV(LV), name, reflMotherLV, isMany, copyNo);		      
    }
		
    return G4PhysicalVolumesPair(pv1, pv2); 			     
  }			   
			     
  //
  //  reflection IS present in transform3D
  //

  if (fVerboseLevel>0)
    G4cout << "scale negative" << G4endl;

  G4VPhysicalVolume* pv1
    = new G4PVPlacement(pureTransform3D, 
                        ReflectLV(LV), name, motherLV, isMany, copyNo);

  G4VPhysicalVolume* pv2 = 0;
  if (G4LogicalVolume* reflMotherLV = GetReflectedLV(motherLV)) {

    // if mother was reflected
    // place the refLV consituent in reflected mother

    pv2 =  new G4PVPlacement(fScale * (pureTransform3D * fScale.inverse()),
                             LV, name, reflMotherLV, isMany, copyNo);
  }			   

  return G4PhysicalVolumesPair(pv1, pv2);  
}					 


//_____________________________________________________________________________
G4PhysicalVolumesPair G4ReflectionFactory::Replicate(const G4String& name, 
	     	                G4LogicalVolume* LV,
			        G4LogicalVolume* motherLV,
                                EAxis axis, 
				G4int nofReplicas, 
		                G4double width,
                                G4double offset)
{
// Creates replica in given mother.
// The result is a pair of physical volumes;
// the second physical volume is a replica in a reflected mother
// - or 0 if mother LV was not reflected.
// ---

  if (fVerboseLevel>0) {
    G4cout << "Replicate " << name << " lv " << LV << " " 
           << LV->GetName() << G4endl;
  }  

  G4VPhysicalVolume* pv1
    = new G4PVReplica(name, LV, motherLV, axis, nofReplicas, width, offset);
 
  G4VPhysicalVolume* pv2 = 0;
  if (G4LogicalVolume* reflMotherLV = GetReflectedLV(motherLV)) {

    // if mother was reflected
    // reflect the LV and replicate it in reflected mother
    
    pv2 = new G4PVReplica(name, ReflectLV(LV), reflMotherLV, 
                          axis, nofReplicas, width, offset); 
  }
		
  return G4PhysicalVolumesPair(pv1, pv2); 			     
}			   
				
			       
//
// private methods
//

//_____________________________________________________________________________
G4LogicalVolume* G4ReflectionFactory::ReflectLV(G4LogicalVolume* LV) 
{
// Gets/creates the reflected solid and logical volume
// and copies + transforms LV daughters.
// ---

  G4LogicalVolume* refLV = GetReflectedLV(LV);

  if (!refLV) {

    // create new (reflected) objects
    refLV = CreateReflectedLV(LV);
			  
    // process daughters  
    ReflectDaughters(LV, refLV);
  }   
  
  return refLV;
}			         

//_____________________________________________________________________________
G4LogicalVolume* G4ReflectionFactory::CreateReflectedLV(G4LogicalVolume* LV) 
{
// Creates the reflected solid and logical volume
// and add the logical volumes pair in the maps.
// ---

  // consistency check
  if (fReflectedLVMap.find(LV) != fReflectedLVMap.end()) {
    G4Exception(
      "G4ReflectionFactory::CreateReflectedLV: called for already reflected volume!");
  }        
				      
  G4VSolid* refSolid 
    = new G4ReflectedSolid(LV->GetSolid()->GetName() + fNameExtension,
                           LV->GetSolid(), fScale);
      
  G4LogicalVolume* refLV
    = new G4LogicalVolume(refSolid, 
                          LV->GetMaterial(),			   	      
			  LV->GetName() + fNameExtension,
			  LV->GetFieldManager(),
			  LV->GetSensitiveDetector(),
			  LV->GetUserLimits());
			  
  fConstituentLVMap[LV] = refLV;
  fReflectedLVMap[refLV] = LV;

  return refLV;			  
}

//_____________________________________________________________________________
void G4ReflectionFactory::ReflectDaughters(G4LogicalVolume* LV, 
                                           G4LogicalVolume* refLV)
{
// Reflects daughters recursively.
// ---

  if (fVerboseLevel>0) {
    G4cout << "G4ReflectionFactory::ReflectDaughters: " 
           << LV->GetNoDaughters() << " of " << LV->GetName() << G4endl;
  }	   

  for (G4int i=0; i<LV->GetNoDaughters(); i++)  {
  
    G4VPhysicalVolume* dPV = LV->GetDaughter(i);
    
    if (! dPV->IsReplicated()) {              
      ReflectPVPlacement(dPV, refLV); 
    }  
    else if (! dPV->GetParameterisation()) {
        ReflectPVReplica(dPV, refLV); 
      }	
      else {                                   
        ReflectPVParameterised(dPV, refLV); 
      }	
  }  
}    
  
//_____________________________________________________________________________
void G4ReflectionFactory::ReflectPVPlacement(G4VPhysicalVolume* dPV, 
                                             G4LogicalVolume* refLV)
{
// Copies and transforms daughter of PVPlacement type of
// a constituent volume into a reflected volume. 
// ---

  G4LogicalVolume* dLV = dPV->GetLogicalVolume();

  // update daughter transformation
  G4Transform3D dt(dPV->GetObjectRotationValue(), dPV->GetObjectTranslation());
  dt = fScale * (dt * fScale.inverse());

  G4LogicalVolume* refDLV;
  
  if (fVerboseLevel>0) 
    G4cout << "Daughter: " << dPV << "  " << dLV->GetName();
  
  if (!IsReflected(dLV)) {

    if (fVerboseLevel>0) 
      G4cout << " will be reflected." << G4endl;

    // create new daughter solid and logical volume
    refDLV = CreateReflectedLV(dLV); 
  
    // create new daughter physical volume
    // with updated transformation

    G4VPhysicalVolume* refDPV
      = new G4PVPlacement(dt, refDLV, dPV->GetName(), refLV, 
	                  dPV->IsMany(), dPV->GetCopyNo()); 
			  
    refLV->AddDaughter(refDPV); 
    
    // recursive call
    ReflectDaughters(dLV, refDLV);   
  } 
  else {
    if (fVerboseLevel>0) 
      G4cout << " will be reconstited." << G4endl;

    refDLV = GetConstituentLV(dLV); 

    G4VPhysicalVolume* refDPV
      = new G4PVPlacement(dt, refDLV, dPV->GetName(), refLV, 
                          dPV->IsMany(), dPV->GetCopyNo()); 
			  
    refLV->AddDaughter(refDPV); 
  }       
}    

//_____________________________________________________________________________
void G4ReflectionFactory::ReflectPVReplica(G4VPhysicalVolume* dPV, 
                                           G4LogicalVolume* refLV)
{
// Copies and transforms daughter of PVReplica type of
// a constituent volume into a reflected volume. 
// ---

  G4LogicalVolume* dLV = dPV->GetLogicalVolume();

  // get replication data
  EAxis axis;
  G4int nofReplicas;
  G4double width;
  G4double offset;
  G4bool consuming;

  dPV->GetReplicationData(axis, nofReplicas, width, offset, consuming);

  G4LogicalVolume* refDLV;
  
  if (fVerboseLevel>0) 
    G4cout << "Daughter: " << dPV << "  " << dLV->GetName();
  
  if (!IsReflected(dLV)) {

    if (fVerboseLevel>0) 
      G4cout << " will be reflected." << G4endl;

    // create new daughter solid and logical volume
    refDLV = CreateReflectedLV(dLV); 
  
    // create new daughter replica

    G4VPhysicalVolume* refDPV
      = new G4PVReplica(dPV->GetName(), refDLV, refLV,
                        axis, nofReplicas, width, offset); 

			  
    refLV->AddDaughter(refDPV); 
    
    // recursive call
    ReflectDaughters(dLV, refDLV);   
  } 
  else {
    if (fVerboseLevel>0) 
      G4cout << " will be reconstited." << G4endl;

    refDLV = GetConstituentLV(dLV); 

    G4VPhysicalVolume* refDPV
      = new G4PVReplica(dPV->GetName(), refDLV, refLV,
                        axis, nofReplicas, width, offset); 
			  
    refLV->AddDaughter(refDPV); 
  }       
}

//_____________________________________________________________________________
void G4ReflectionFactory::ReflectPVParameterised(G4VPhysicalVolume* dPV, 
                                           G4LogicalVolume* refLV)
{
// Not implemented.
// Should copy and transform daughter of PVReplica type of
// a constituent volume into a reflected volume. 
// ---

  G4Exception(
    "G4ReflectionFactory: Parameterised volumes cannot be reflected yet.");

}

//_____________________________________________________________________________
G4LogicalVolume* G4ReflectionFactory::GetConstituentLV(
                                         G4LogicalVolume* reflLV) const
{				      
// Returns the consituent volume of the given reflected volume,
// 0 if the given reflected volume was not found.
// ---

  LogicalVolumesMapIterator it = fReflectedLVMap.find(reflLV);

  if (it == fReflectedLVMap.end()) return 0;

  return (*it).second;
}			  

//_____________________________________________________________________________
G4LogicalVolume* G4ReflectionFactory::GetReflectedLV(
                                         G4LogicalVolume* lv) const
{				      
// Returns the reflected volume of the given consituent volume,
// 0 if the given volume was not reflected.
// ---

  LogicalVolumesMapIterator it = fConstituentLVMap.find(lv);

  if (it == fConstituentLVMap.end()) return 0;

  return (*it).second;
}			  

//_____________________________________________________________________________
G4bool G4ReflectionFactory::IsConstituent(G4LogicalVolume* lv) const
{
// Returns true if the given volume has been already reflected
// (is in the map of constituent volumes).
// ---

  return (fConstituentLVMap.find(lv) != fConstituentLVMap.end());
}  

//_____________________________________________________________________________
G4bool G4ReflectionFactory::IsReflected(G4LogicalVolume* lv) const
{
// Returns true if the given volume is a reflected volume
// (is in the map reflected  volumes).
// ---

  return (fReflectedLVMap.find(lv) != fReflectedLVMap.end());
}  

//_____________________________________________________________________________
G4bool G4ReflectionFactory::IsReflection(const G4Scale3D& scale) const
{
// Returns true if the scale is negative, false otherwise.
// ---

  if (scale(0,0)*scale(1,1)*scale(2,2) < 0.)
    return true;
  else 
    return false;  
}

//_____________________________________________________________________________
void G4ReflectionFactory::PrintConstituentLVMap()
{
// temporary - for debugging purpose
// ---

  LogicalVolumesMapIterator it;
  for (it = fConstituentLVMap.begin(); it != fConstituentLVMap.end(); it++) {
    G4cout << "lv: " << (*it).first << "  lv_refl: " << (*it).second << G4endl;	  
  }
  G4cout << G4endl;
}  

//_____________________________________________________________________________
void G4ReflectionFactory::CheckScale(const G4Scale3D& scale) const
{
// Check if scale correspond to fScale,
// if not give exception.
// ---

  if (!IsReflection(scale)) return;
  
  G4double diff = 0.;
  for (G4int i=0; i<4; i++)
    for (G4int j=0; j<4; j++) 
      diff += abs(scale(i,j) - fScale(i,j));  

  if (diff > kCarTolerance)
    G4Exception("G4ReflectionFactory::CheckScale: unexpected scale has occured.");
}    

void  G4ReflectionFactory::SetVerboseLevel(G4int verboseLevel)
{
  fVerboseLevel = verboseLevel;
}
				  
G4int G4ReflectionFactory::GetVerboseLevel() const
{
  return fVerboseLevel;
}

/*
  // placement with decomposed transformation

  G4VPhysicalVolume* pv1
    =  new G4PVPlacement(new G4RotationMatrix(rotation.getRotation().inverse()),
                         translation.getTranslation(),
    	 		 refLV, name, motherLV, isMany, copyNo);
*/			 
