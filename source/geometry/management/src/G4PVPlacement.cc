// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVPlacement.cc,v 1.3 2000-11-20 17:31:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVPlacement Implementation

#include "G4PVPlacement.hh"
// #include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"    

G4PVPlacement::G4PVPlacement(G4RotationMatrix *pRot,
			     const G4ThreeVector &tlate,
			     const G4String& pName,
			     G4LogicalVolume *pLogical,
			     G4VPhysicalVolume *pMother,
			     G4bool pMany,
			     G4int pCopyNo) :
  G4VPhysicalVolume(pRot,tlate,pName,pLogical,pMother),
  fmany(pMany), fallocatedRotM(false), fcopyNo(pCopyNo)
{
}

G4PVPlacement::G4PVPlacement(const G4Transform3D &Transform3D,
			     const G4String& pName,
			     G4LogicalVolume *pLogical,
			     G4VPhysicalVolume *pMother,
			     G4bool pMany,
			     G4int pCopyNo) :
  G4VPhysicalVolume( NewPtrRotMatrix(Transform3D.getRotation().inverse()),
		    Transform3D.getTranslation(),
		    pName,pLogical,pMother),
  fmany(pMany), fcopyNo(pCopyNo)
{
  fallocatedRotM= (this->GetRotation() != 0);
}

//
//   The logical volume of the mother is utilised (not the physical)
//   [ The physical volume is needed, known and set only at tracking time. ]
//

G4PVPlacement::G4PVPlacement(G4RotationMatrix *pRot,
			     const G4ThreeVector &tlate,
			     G4LogicalVolume *pCurrentLogical,
			     const G4String& pName,
			     G4LogicalVolume *pMotherLogical,
			     G4bool pMany,
			     G4int pCopyNo) :
  G4VPhysicalVolume(pRot,tlate,pName,pCurrentLogical,0),
  fmany(pMany), fallocatedRotM(false), fcopyNo(pCopyNo)
{
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
}


G4PVPlacement::G4PVPlacement( const G4Transform3D &Transform3D,
			      G4LogicalVolume *pCurrentLogical,
			      const G4String& pName,
			      G4LogicalVolume *pMotherLogical,
			      G4bool pMany,
			      G4int pCopyNo):
    G4VPhysicalVolume( 0,
		      Transform3D.getTranslation(),
		      pName,
                      pCurrentLogical,
                      0),
    fmany(pMany), fcopyNo(pCopyNo)
{
  this->SetRotation( NewPtrRotMatrix(Transform3D.getRotation().inverse()) );
  fallocatedRotM= (this->GetRotation() != 0);
  
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
}



G4PVPlacement::~G4PVPlacement()
{
  if( fallocatedRotM ){ delete frot; }
}

G4bool G4PVPlacement::IsMany() const
{
    return fmany; 
}

G4int G4PVPlacement::GetCopyNo() const
{
    return fcopyNo;
}


void  G4PVPlacement::SetCopyNo(G4int newCopyNo)
{
    fcopyNo= newCopyNo;
}


G4bool G4PVPlacement::IsReplicated() const
{
    return false;
}

G4VPVParameterisation* G4PVPlacement::GetParameterisation() const
{
    return 0;
}

void G4PVPlacement::GetReplicationData(EAxis& axis,
                                   G4int& nReplicas,
				   G4double& width,
                                   G4double& offset,
                                   G4bool& consuming) const
{
// No-operations
}

void G4PVPlacement::Setup(G4VPhysicalVolume *pMother)
{
    SetMother(pMother);
}

//  
// Auxiliary function for 2nd & 4th constructors (the ones with G4Transform3D)
//  Creates a new RotMatrix on the heap (using "new") and copies 
//  its argument into it.
//
//  No entity is currently responsible to delete this memory. 
//  This is a memory leak.  <-- FIXME
//

G4RotationMatrix* G4PVPlacement::NewPtrRotMatrix(const G4RotationMatrix &RotMat)
{
    G4RotationMatrix *pRotMatrix; 
    if ( RotMat.isIdentity() )
       pRotMatrix= 0;
    else{
       pRotMatrix= new G4RotationMatrix(RotMat);
    }
    // fallocatedRotM= ! (RotMat.isIdentity());
    
    return pRotMatrix;
}
