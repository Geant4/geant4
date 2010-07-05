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
// $Id: G4PVReplica.cc,v 1.8 2010-07-05 13:29:12 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVReplica Implementation
//
// ----------------------------------------------------------------------

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
    fcopyNo(-1), fRegularVolsId(0)
{
  if ((!pMother) || (!pMother->GetLogicalVolume()))
  {
    G4cerr << "ERROR - NULL pointer specified as mother volume." << G4endl
           << "        The world volume cannot be sliced or parameterised !"
           << G4endl;
    G4Exception("G4PVReplica::G4PVReplica()", "InvalidSetup", FatalException,
                "NULL pointer to mother. World volume cannot be sliced !");
    return;
  }
  G4LogicalVolume* motherLogical = pMother->GetLogicalVolume();
  if (pLogical == motherLogical)
  {
    G4Exception("G4PVReplica::G4PVReplica()", "InvalidSetup",
                FatalException, "Cannot place a volume inside itself!");
    return;
  }
  SetMotherLogical(motherLogical);
  motherLogical->AddDaughter(this);
  if (motherLogical->GetNoDaughters() != 1)
  {
    G4cerr << "ERROR - A replica or parameterised volume must be" << G4endl
           << "        the only daughter of a given mother volume !" << G4endl
           << "           Mother physical volume: " << pMother->GetName() << G4endl
           << "           Replicated volume: " << pName << G4endl;
    G4Exception("G4PVReplica::G4PVReplica()", "InvalidSetup", FatalException,
                "Replica or parameterised volume must be the only daughter !");
    return;
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
    fcopyNo(-1), fRegularVolsId(0)
{
  if (!pMotherLogical)
  {
    G4cerr << "ERROR - NULL pointer specified as mother volume for "
           << pName << "." << G4endl;
    G4Exception("G4PVReplica::G4PVReplica()", "NullPointer", FatalException,
                "Invalid setup. NULL pointer specified as mother !");
    return;
  }
  if (pLogical == pMotherLogical)
  {
    G4Exception("G4PVReplica::G4PVReplica()", "InvalidSetup",
                FatalException, "Cannot place a volume inside itself!");
    return;
  }
  pMotherLogical->AddDaughter(this);
  SetMotherLogical(pMotherLogical);
  if (pMotherLogical->GetNoDaughters() != 1)
  {
    G4cerr << "ERROR - A replica or parameterised volume must be" << G4endl
           << "        the only daughter of a given mother volume !" << G4endl
           << "           Mother logical volume: " << pMotherLogical->GetName() << G4endl
           << "           Replicated volume: " << pName << G4endl;
    G4Exception("G4PVReplica::G4PVReplica()", "InvalidSetup", FatalException,
                "Replica or parameterised volume must be the only daughter !");
    return;
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
    G4Exception("G4PVReplica::CheckAndSetParameters()", "WrongArgumentValue",
                FatalException, "Illegal number of replicas.");
  }
  fnReplicas=nReplicas;
  if (width<0)
  {
    G4Exception("G4PVReplica::CheckAndSetParameters()", "WrongArgumentValue",
                FatalException, "Width must be positive.");
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
        G4Exception("G4PVReplica::CheckAndSetParameters()", "FatalError",
                    FatalException, "Rotation matrix allocation failed.");
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
      G4Exception("G4PVReplica::CheckAndSetParameters()", "WrongArgumentValue",
                  FatalException, "Unknown axis of replication.");
      break;
  }
}

G4PVReplica::G4PVReplica( __void__& a )
  : G4VPhysicalVolume(a), faxis(kZAxis), fnReplicas(0), fwidth(0.),
    foffset(0.), fcopyNo(-1), fRegularStructureCode(0), fRegularVolsId(0)
{
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

G4bool G4PVReplica::IsRegularStructure() const
{
  return (fRegularVolsId!=0); 
}

G4int  G4PVReplica::GetRegularStructureId() const
{
  return fRegularVolsId; 
}

void   G4PVReplica::SetRegularStructureId( G4int Code )
{
  fRegularVolsId= Code; 
} 
