// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVReplica.cc,v 1.2 1999-12-15 14:49:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVReplica Implementation

#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"

G4PVReplica::G4PVReplica(const G4String& pName,
			 G4LogicalVolume* pLogical,
			 G4VPhysicalVolume* pMother,
                         const EAxis pAxis,
                         const G4int nReplicas,
		         const G4double width,
                         const G4double offset) : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,pMother), fcopyNo(-1)
{
  CheckAndSetParameters (pAxis, nReplicas, width, offset);
}

G4PVReplica::G4PVReplica(const G4String& pName,
			 G4LogicalVolume* pLogical,
			 G4LogicalVolume* pMotherLogical,
                         const EAxis pAxis,
                         const G4int nReplicas,
		         const G4double width,
                         const G4double offset) : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
  CheckAndSetParameters (pAxis, nReplicas, width, offset);
}

void G4PVReplica::CheckAndSetParameters (const EAxis pAxis,
				    const G4int nReplicas,
				    const G4double width,
				    const G4double offset) {
    
    if (nReplicas<1)
	{
        G4Exception("G4PVReplica::G4PVReplica illegal no of replicas");
	}
    fnReplicas=nReplicas;
    if (width<0)
	{
        G4Exception("G4PVReplica::G4PVReplica Width must be positive");
	}
    fwidth=width;

    foffset=offset;

    faxis=pAxis;

// Create rotation matrix for phi axis case & check axis is valid
    G4RotationMatrix *pRMat;
    switch (faxis)
	{
	case kPhi:
	  pRMat=new G4RotationMatrix();
	  if (!pRMat)
	    {
	      G4Exception("G4PVReplica::G4PVReplica rotation matrix alloc failed");
	    }
	  SetRotation(pRMat);
	  break;
        case kRho:
	case kXAxis:
	case kYAxis:
	case kZAxis:
	  break;
	default:
	  G4Exception("G4PVReplica::G4PVReplica Unknown axis");
	  break;
	}

}

G4PVReplica::~G4PVReplica()
{
     if (faxis==kPhi)
 	{
 	delete GetRotation();
 	}
}

G4bool G4PVReplica::IsMany() const
{
    return false; 
}

G4int G4PVReplica::GetCopyNo() const
{
    // return 0;
    return fcopyNo;
}

void  G4PVReplica::SetCopyNo(G4int newCopyNo)
{
    fcopyNo= newCopyNo;
}

G4bool G4PVReplica::IsReplicated() const
{
    return true;
}

G4VPVParameterisation* G4PVReplica::GetParameterisation() const
{
    return 0;
}

void G4PVReplica::GetReplicationData(EAxis& axis,
                                   G4int& nReplicas,
				   G4double& width,
                                   G4double& offset,
                                   G4bool& consuming) const
{
    axis=faxis;
    nReplicas=fnReplicas;
    width=fwidth;
    offset=foffset;
    consuming=true;
}

void G4PVReplica::Setup(G4VPhysicalVolume *pMother)
{
    SetMother(pMother);
}




