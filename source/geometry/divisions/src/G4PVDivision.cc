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
// $Id: G4PVDivision.cc,v 1.1 2003-06-16 15:11:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4PVDivision Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "G4PVDivision.hh"
#include "G4LogicalVolume.hh"
#include "G4ParameterisationBox.hh"
#include "G4ParameterisationTrd.hh"
#include "G4ParameterisationPolycone.hh"
#include "G4ParameterisationCons.hh"
#include "G4VSolid.hh"

//--------------------------------------------------------------------------
G4PVDivision::G4PVDivision(const G4String& pName,
                                 G4LogicalVolume* pLogical,
                                 G4VPhysicalVolume* pMother,
                           const EAxis pAxis,
                           const G4int nReplicas,
                           const G4double width,
                           const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,pMother),
    fcopyNo(-1)
{
  SetParameterisation( pMother->GetLogicalVolume(), pAxis,
                       nReplicas, width, offset, DivNDIVandWIDTH );
  CheckAndSetParameters (pAxis, nReplicas, width, offset, DivNDIVandWIDTH);
}

//--------------------------------------------------------------------------
G4PVDivision::G4PVDivision(const G4String& pName,
                                 G4LogicalVolume* pLogical,
                                 G4LogicalVolume* pMotherLogical,
                           const EAxis pAxis,
                           const G4int nReplicas,
                           const G4double width,
                           const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0),
    fcopyNo(-1)
{
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
  SetParameterisation( pMotherLogical, pAxis, nReplicas,
                       width, offset, DivNDIVandWIDTH );
  CheckAndSetParameters (pAxis, nReplicas, width, offset, DivNDIVandWIDTH);
}

//--------------------------------------------------------------------------
G4PVDivision::G4PVDivision(const G4String& pName,
                                 G4LogicalVolume* pLogical,
                                 G4LogicalVolume* pMotherLogical,
                           const EAxis pAxis,
                           const G4int nReplicas,
                           const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0),
    fcopyNo(-1)
{
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
  SetParameterisation( pMotherLogical, pAxis, nReplicas, 0., offset, DivNDIV );
  CheckAndSetParameters (pAxis, nReplicas, 0., offset, DivNDIV );
}

//--------------------------------------------------------------------------
G4PVDivision::G4PVDivision(const G4String& pName,
                                 G4LogicalVolume* pLogical,
                                 G4LogicalVolume* pMotherLogical,
                           const EAxis pAxis,
                           const G4double width,
                           const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0),
    fcopyNo(-1)
{
  if (pMotherLogical) pMotherLogical->AddDaughter(this);
  SetParameterisation( pMotherLogical, pAxis, 0, width, offset, DivWIDTH );
  CheckAndSetParameters (pAxis, 0, width, offset, DivWIDTH );
}

//--------------------------------------------------------------------------
void G4PVDivision::CheckAndSetParameters (const EAxis pAxis,
                                          const G4int nReplicas,
                                          const G4double width,
                                          const G4double offset,
                                                DivisionType divType ) 
{
  if( divType == DivWIDTH )
  {
    fnReplicas = fparam->GetNoDiv();
  }
  else
  {
    fnReplicas = nReplicas;
  }
  if (fnReplicas < 1 )
  {
    G4Exception("G4PVDivision::G4PVDivision() - Illegal number of replicas!");
  }

  if( divType != DivNDIV)
  {
    fwidth = fparam->GetWidth();
  }
  else
  {
    fwidth = width;
  }
  if( fwidth < 0 )
  {
    G4Exception("G4PVDivision::G4PVDivision() - Width must be positive!");
  }
  
  foffset=offset;
  
  //!!!!! axis has to be x/y/z in G4VoxelLimits::GetMinExtent
  //
  if( pAxis == kRho || pAxis == kRadial3D || pAxis == kPhi )
  {
    faxis = kZAxis;
  }
  else
  {
    faxis = pAxis;
  }
  
  // Create rotation matrix: for phi axis it will be changed
  // in G4VPVParameterisation::ComputeTransformation, for others
  // it will stay the unity
  //
  G4RotationMatrix *pRMat = new G4RotationMatrix();
  SetRotation(pRMat);
  
  switch (faxis)
  {
    case kPhi:
      break;
    case kRho:
    case kXAxis:
    case kYAxis:
    case kZAxis:
      break;
    default:
      G4Exception("G4PVDivision::G4PVDivision() - Unknown axis.");
      break;
  }
}

//--------------------------------------------------------------------------
G4PVDivision::~G4PVDivision()
{
  delete GetRotation();
}

//--------------------------------------------------------------------------
G4bool G4PVDivision::IsParameterised() const
{ 
  return true;
}

//--------------------------------------------------------------------------
G4bool G4PVDivision::IsMany() const
{
  return false; 
}

//--------------------------------------------------------------------------
G4int G4PVDivision::GetCopyNo() const
{
  return fcopyNo;
}

//--------------------------------------------------------------------------
void  G4PVDivision::SetCopyNo(G4int newCopyNo)
{
  fcopyNo= newCopyNo;
}

//--------------------------------------------------------------------------
G4bool G4PVDivision::IsReplicated() const
{
  return true;
}

//--------------------------------------------------------------------------
G4VPVParameterisation* G4PVDivision::GetParameterisation() const
{
  return fparam;
}

//--------------------------------------------------------------------------
void G4PVDivision::GetReplicationData(EAxis& axis,
                                      G4int& nReplicas,
                                      G4double& width,
                                      G4double& offset,
                                      G4bool& consuming ) const
{
  axis=faxis;
  nReplicas=fnReplicas;
  width=fwidth;
  offset=foffset;
  consuming=false;
}

//--------------------------------------------------------------------------
void G4PVDivision::Setup(G4VPhysicalVolume *pMother)
{
  SetMother(pMother);
}


//--------------------------------------------------------------------------
void G4PVDivision::SetParameterisation( G4LogicalVolume* motherLogical,
                                  const EAxis axis,
                                  const G4int nReplicas,
                                  const G4double width,
                                  const G4double offset,
                                        DivisionType divType )
{
  // Check that solid is compatible with mother solid and axis of division   
  // CheckSolid( solid, motherSolid );
  // G4cout << " Axis " << axis << G4endl;

  G4VSolid* mSolid = motherLogical->GetSolid();
  G4String mSolidType = mSolid->GetEntityType();

  // Parameterisation type depend of mother solid type and axis of division
  //
  if( mSolidType == "G4Box" )
  {
    switch( axis )
    {
      case kXAxis:
        fparam = new G4ParameterisationBoxX( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      case kYAxis:
        fparam = new G4ParameterisationBoxY( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationBoxZ( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
        break;
    }
  }
  else if( mSolidType == "G4Trd" )
  { 
    switch( axis )
    {
      case kXAxis:
        fparam = new G4ParameterisationTrdX( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      case kYAxis:
        fparam = new G4ParameterisationTrdY( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationTrdZ( axis, nReplicas, width,
                                             offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
        break;
    }
  }
//  else if( mSolidType == "G4Trap" )
//  {
//  }
  else if( mSolidType == "G4Polycone" )
  {
    switch( axis )
    {
      case kRho:
        fparam = new G4ParameterisationPolyconeRho( axis, nReplicas, width,
                                                    offset, mSolid, divType );
        break;
      case kPhi:
        fparam = new G4ParameterisationPolyconePhi( axis, nReplicas, width,
                                                    offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationPolyconeZ( axis, nReplicas, width,
                                                  offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
      break;
    }
  }
  else if( mSolidType == "G4Cons" )
  {
    switch( axis )
    {
      case kRho:
        fparam = new G4ParameterisationConsRho( axis, nReplicas, width,
                                                offset, mSolid, divType );
        break;
      case kPhi:
        fparam = new G4ParameterisationConsPhi( axis, nReplicas, width,
                                                offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationConsZ( axis, nReplicas, width,
                                              offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
        break;
    }
  }
  else
  {
    G4cerr << "ERROR - G4PVDivision::SetParameterisation()" << G4endl
           << "        Divisions for" << mSolidType
           << " not implemented." << G4endl;
    G4Exception("G4PVDivision - Solid type not supported.");
  }
}

//--------------------------------------------------------------------------
void G4PVDivision::ErrorInAxis( EAxis axis, G4VSolid* solid )
{
  G4String error = "G4PVDivision - Error in axis: trying to divide solid "
                 + solid->GetName() + " of type " + solid->GetEntityType()
                 + " along axis ";
  switch( axis )
  {
    case kXAxis:
      error += "X.";
      break;
    case kYAxis:
      error += "Y.";
      break;
    case kZAxis:
      error += "Z.";
      break;
    case kRho:
      error += "Rho.";
      break;
    case kRadial3D:
      error += "Radial3D.";
      break;
    case kPhi:
      error += "Phi.";
      break;
    default:
      break;
  }
  G4Exception( error );
}
