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
// $Id: G4PVReplica.cc,v 1.9 2002-10-23 16:26:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVReplica Implementation
//
// ********************************************************************

#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"

G4PVReplica::G4PVReplica( const G4String& pName,
                                G4LogicalVolume* pLogical,
                                G4VPhysicalVolume* pMother,
                          const EAxis pAxis,
                          const G4int nReplicas,
                          const G4double width,
                          const G4double offset )
  : G4VPhysicalVolume(0, G4ThreeVector(), pName, pLogical, pMother),
    fcopyNo(-1)
{
  if ((!pMother) || (!pMother->GetLogicalVolume()))
  {
    G4cout << "ERROR - NULL pointer specified as mother volume." << G4endl
           << "        The world volume cannot be sliced or parameterised !"
           << G4endl;
    G4Exception("G4PVReplica::G4PVReplica() - Not a valid geometry setup !");
  }
  if (pMother->GetLogicalVolume()->GetNoDaughters() != 1)
  {
    G4cout << "ERROR - A replica or parameterised volume must be" << G4endl
           << "        the only daughter of a given mother volume !" << G4endl
           << "           Mother physical volume: " << pMother->GetName() << G4endl
           << "           Replicated volume: " << pName << G4endl;
    G4Exception("G4PVReplica::G4PVReplica() - Not a valid geometry setup !");
  }
  CheckAndSetParameters (pAxis, nReplicas, width, offset);
}

G4PVReplica::G4PVReplica( const G4String& pName,
                                G4LogicalVolume* pLogical,
                                G4LogicalVolume* pMotherLogical,
                          const EAxis pAxis,
                          const G4int nReplicas,
                          const G4double width,
                          const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0),
    fcopyNo(-1)
{
  if (!pMotherLogical)
  {
    G4cout << "ERROR - NULL pointer specified as mother volume for "
           << pName << "." << G4endl;
    G4Exception("G4PVReplica::G4PVReplica() - Not a valid geometry setup !");
  }
  pMotherLogical->AddDaughter(this);
  if (pMotherLogical->GetNoDaughters() != 1)
  {
    G4cout << "ERROR - A replica or parameterised volume must be" << G4endl
           << "        the only daughter of a given mother volume !" << G4endl
           << "           Mother logical volume: " << pMotherLogical->GetName() << G4endl
           << "           Replicated volume: " << pName << G4endl;
    G4Exception("G4PVReplica::G4PVReplica() - Not a valid geometry setup");
  }
  CheckAndSetParameters (pAxis, nReplicas, width, offset);
}

void G4PVReplica::CheckAndSetParameters( const EAxis pAxis,
                                         const G4int nReplicas,
                                         const G4double width,
                                         const G4double offset)
{    
  if (nReplicas<1)
  {
    G4Exception("G4PVReplica::CheckAndSetParameters() - Illegal number of replicas");
  }
  fnReplicas=nReplicas;
  if (width<0)
  {
    G4Exception("G4PVReplica::CheckAndSetParameters() - Width must be positive");
  }
  fwidth  = width;
  foffset = offset;
  faxis   = pAxis;

  // Create rotation matrix for phi axis case & check axis is valid
  //
  G4RotationMatrix* pRMat=0;
  switch (faxis)
  {
    case kPhi:
      pRMat=new G4RotationMatrix();
      if (!pRMat)
      {
        G4Exception("G4PVReplica::CheckAndSetParameters() - Rotation matrix alloc failed");
      }
      SetRotation(pRMat);
      break;
    case kRho:
    case kXAxis:
    case kYAxis:
    case kZAxis:
    case kUndefined:
      break;
    default:
      G4Exception("G4PVReplica::CheckAndSetParameters() - Unknown axis");
      break;
  }
}

G4PVReplica::~G4PVReplica()
{
  if ( faxis==kPhi )
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
  return fcopyNo;
}

void  G4PVReplica::SetCopyNo(G4int newCopyNo)
{
  fcopyNo = newCopyNo;
}

G4bool G4PVReplica::IsReplicated() const
{
  return true;
}

G4bool G4PVReplica::IsParameterised() const
{
  return false;
}

G4VPVParameterisation* G4PVReplica::GetParameterisation() const
{
  return 0;
}

G4int G4PVReplica::GetMultiplicity() const
{
  return fnReplicas;
}

void G4PVReplica::GetReplicationData( EAxis& axis,
                                      G4int& nReplicas,
                                      G4double& width,
                                      G4double& offset,
                                      G4bool& consuming ) const
{
  axis = faxis;
  nReplicas = fnReplicas;
  width = fwidth;
  offset = foffset;
  consuming = true;
}

void G4PVReplica::Setup(G4VPhysicalVolume *pMother)
{
  SetMother(pMother);
}
