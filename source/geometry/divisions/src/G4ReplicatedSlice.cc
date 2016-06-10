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
// $Id: G4ReplicatedSlice.cc 89815 2015-04-30 14:50:22Z gcosmo $
//
// --------------------------------------------------------------------

#include "G4ReplicatedSlice.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ReflectedSolid.hh"
#include "G4ParameterisationBox.hh"
#include "G4ParameterisationTubs.hh"
#include "G4ParameterisationCons.hh"
#include "G4ParameterisationTrd.hh"
#include "G4ParameterisationPara.hh"
#include "G4ParameterisationPolycone.hh"
#include "G4ParameterisationPolyhedra.hh"

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4LogicalVolume* pMotherLogical,
                                     const EAxis pAxis,
                                     const G4int nDivs,
                                     const G4double width,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, nDivs, width, half_gap, offset,
                        DivNDIVandWIDTH, pMotherLogical, pLogical);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4LogicalVolume* pMotherLogical,
                                     const EAxis pAxis,
                                     const G4int nDivs,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, nDivs, 0., half_gap, offset,
                        DivNDIV, pMotherLogical, pLogical);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4LogicalVolume* pMotherLogical,
                                     const EAxis pAxis,
                                     const G4double width,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, 0, width, half_gap, offset,
                        DivWIDTH, pMotherLogical, pLogical);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4VPhysicalVolume* pMotherPhysical,
                                     const EAxis pAxis,
                                     const G4int nDivs,
                                     const G4double width,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, nDivs, width, half_gap, offset,
      DivNDIVandWIDTH, pMotherPhysical->GetLogicalVolume(), pLogical);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4VPhysicalVolume* pMotherPhysical,
                                     const EAxis pAxis,
                                     const G4int nDivs,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, nDivs, 0., half_gap, offset,
      DivNDIV, pMotherPhysical->GetLogicalVolume(), pLogical);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::G4ReplicatedSlice(const G4String& pName,
                                           G4LogicalVolume* pLogical,
                                           G4VPhysicalVolume* pMotherPhysical,
                                     const EAxis pAxis,
                                     const G4double width,
                                     const G4double half_gap,
                                     const G4double offset )
  : G4VPhysicalVolume(0,G4ThreeVector(),pName,pLogical,0), fcopyNo(-1)
{
  CheckAndSetParameters(pAxis, 0, width, half_gap, offset,
      DivWIDTH, pMotherPhysical->GetLogicalVolume(), pLogical);
}

//--------------------------------------------------------------------------
void
G4ReplicatedSlice::CheckAndSetParameters( const EAxis pAxis,
                                          const G4int nDivs,
                                          const G4double width,
                                          const G4double half_gap,
                                          const G4double offset, 
                                                DivisionType divType,
                                                G4LogicalVolume* pMotherLogical,
                                          const G4LogicalVolume* pLogical )
{
  if(!pMotherLogical)
  {
    std::ostringstream message;
    message << "Invalid setup." << G4endl
            << "NULL pointer specified as mother! Volume: " << GetName();
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                FatalException, message);
  }
  if(pLogical == pMotherLogical)
  {
    std::ostringstream message;
    message << "Invalid setup." << G4endl
            << "Cannot place a volume inside itself! Volume: " << GetName();
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                FatalException, message);
  }

  //----- Check that mother solid is of the same type as
  //      daughter solid (otherwise, the corresponding
  //      Parameterisation::ComputeDimension() will not be called)
  //
  G4String msolType = pMotherLogical->GetSolid()->GetEntityType();
  G4String dsolType = pLogical->GetSolid()->GetEntityType();
  if( msolType != dsolType && ( msolType != "G4Trd" || dsolType != "G4Trap" ) )
  {
    std::ostringstream message;
    message << "Invalid setup." << G4endl
            << "Incorrect solid type for division of volume: "
            << GetName() << G4endl
            << "    It is: " << msolType
            << ", while it should be: " << dsolType;
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()",
                "GeomDiv0002", FatalException, message);
  }

  pMotherLogical->AddDaughter(this);
  SetMotherLogical(pMotherLogical);
  SetParameterisation(pMotherLogical, pAxis, nDivs,
                      width, half_gap, offset, divType);

  if( divType == DivWIDTH )
  {
    fnReplicas = fparam->GetNoDiv();
  }
  else
  {
    fnReplicas = nDivs;
  }
  if (fnReplicas < 1 )
  {
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                FatalException, "Illegal number of replicas!");
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
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                FatalException, "Width must be positive!");
  }
  if( fwidth < 2.*half_gap )
  {
    G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                FatalException, "Half_gap is too large!");
  }
  
  foffset = offset;
  fdivAxis = pAxis;

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
  
  switch (faxis)
  {
    case kPhi:
    case kRho:
    case kXAxis:
    case kYAxis:
    case kZAxis:
      break;
    default:
      G4Exception("G4ReplicatedSlice::CheckAndSetParameters()", "GeomDiv0002",
                  FatalException, "Unknown axis of replication.");
      break;
  }

  // Create rotation matrix: for phi axis it will be changed
  // in G4VPVParameterisation::ComputeTransformation, for others
  // it will stay the unity
  //
  G4RotationMatrix *pRMat = new G4RotationMatrix();
  SetRotation(pRMat);
}

//--------------------------------------------------------------------------
G4ReplicatedSlice::~G4ReplicatedSlice()
{
  delete GetRotation();
}

//--------------------------------------------------------------------------
EAxis G4ReplicatedSlice::GetDivisionAxis() const
{
  return fdivAxis;
}

//--------------------------------------------------------------------------
G4bool G4ReplicatedSlice::IsParameterised() const
{ 
  return true;
}

//--------------------------------------------------------------------------
G4bool G4ReplicatedSlice::IsMany() const
{
  return false; 
}

//--------------------------------------------------------------------------
G4int G4ReplicatedSlice::GetCopyNo() const
{
  return fcopyNo;
}

//--------------------------------------------------------------------------
void  G4ReplicatedSlice::SetCopyNo(G4int newCopyNo)
{
  fcopyNo= newCopyNo;
}

//--------------------------------------------------------------------------
G4bool G4ReplicatedSlice::IsReplicated() const
{
  return true;
}

//--------------------------------------------------------------------------
G4VPVParameterisation* G4ReplicatedSlice::GetParameterisation() const
{
  return fparam;
}

//--------------------------------------------------------------------------
void G4ReplicatedSlice::GetReplicationData(EAxis& axis,
                                           G4int& nDivs,
                                           G4double& width,
                                           G4double& offset,
                                           G4bool& consuming ) const
{
  axis=faxis;
  nDivs=fnReplicas;
  width=fwidth;
  offset=foffset;
  consuming=false;
}


//--------------------------------------------------------------------------
void G4ReplicatedSlice::SetParameterisation( G4LogicalVolume* motherLogical,
                                       const EAxis axis,
                                       const G4int nDivs,
                                       const G4double width,
                                       const G4double half_gap,
                                       const G4double offset,
                                             DivisionType divType )
{
  G4VSolid* mSolid = motherLogical->GetSolid();
  G4String mSolidType = mSolid->GetEntityType();
  fparam = 0;

  // If the solid is a reflected one, update type to its
  // real constituent solid.
  // 
  if (mSolidType == "G4ReflectedSolid")
  {
      mSolidType = ((G4ReflectedSolid*)mSolid)->GetConstituentMovedSolid()
                 ->GetEntityType(); 
  }    

  // Parameterisation type depend of mother solid type and axis of division
  //
  if( mSolidType == "G4Box" )
  {
    switch( axis )
    {
      case kXAxis:
        fparam = new G4ParameterisationBoxX( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kYAxis:
        fparam = new G4ParameterisationBoxY( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationBoxZ( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
        break;
    }
  }
  else if( mSolidType == "G4Tubs" )
  {
    switch( axis )
    {
      case kRho:
        fparam = new G4ParameterisationTubsRho( axis, nDivs, width,
                                                offset, mSolid, divType );
        break;
      case kPhi:
        fparam = new G4ParameterisationTubsPhi( axis, nDivs, width,
                                                offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationTubsZ( axis, nDivs, width,
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
        fparam = new G4ParameterisationConsRho( axis, nDivs, width,
                                                offset, mSolid, divType );
        break;
      case kPhi:
        fparam = new G4ParameterisationConsPhi( axis, nDivs, width,
                                                offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationConsZ( axis, nDivs, width,
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
        fparam = new G4ParameterisationTrdX( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kYAxis:
        fparam = new G4ParameterisationTrdY( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationTrdZ( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      default:
        ErrorInAxis( axis, mSolid );
        break;
    }
  }
  else if( mSolidType == "G4Para" )
  { 
    switch( axis )
    {
      case kXAxis:
        fparam = new G4ParameterisationParaX( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kYAxis:
        fparam = new G4ParameterisationParaY( axis, nDivs, width,
                                             offset, mSolid, divType );
        break;
      case kZAxis:
        fparam = new G4ParameterisationParaZ( axis, nDivs, width,
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
//  else if( mSolidType == "G4Polycone" )
//  {
//    switch( axis )
//    {
//      case kRho:
//        fparam = new G4ParameterisationPolyconeRho( axis, nDivs, width,
//                                                    offset, mSolid, divType );
//        break;
//      case kPhi:
//        fparam = new G4ParameterisationPolyconePhi( axis, nDivs, width,
//                                                    offset, mSolid, divType );
//        break;
//      case kZAxis:
//        fparam = new G4ParameterisationPolyconeZ( axis, nDivs, width,
//                                                  offset, mSolid, divType );
//        break;
//      default:
//        ErrorInAxis( axis, mSolid );
//      break;
//    }
//  }
//  else if( mSolidType == "G4Polyhedra" )
//  {
//    switch( axis )
//    {
//      case kRho:
//        fparam = new G4ParameterisationPolyhedraRho( axis, nDivs, width,
//                                                    offset, mSolid, divType );
//        break;
//      case kPhi:
//        fparam = new G4ParameterisationPolyhedraPhi( axis, nDivs, width,
//                                                    offset, mSolid, divType );
//        break;
//      case kZAxis:
//        fparam = new G4ParameterisationPolyhedraZ( axis, nDivs, width,
//                                                  offset, mSolid, divType );
//        break;
//      default:
//        ErrorInAxis( axis, mSolid );
//      break;
//    }
//  }
  else
  {
    std::ostringstream message;
    message << "Solid type not supported: " << mSolidType << "." << G4endl
            << "Divisions for " << mSolidType << " not implemented.";
    G4Exception("G4ReplicatedSlice::SetParameterisation()", "GeomDiv0001",
                FatalException, message);
  }

  fparam->SetHalfGap(half_gap);
}

//--------------------------------------------------------------------------
void G4ReplicatedSlice::ErrorInAxis( EAxis axis, G4VSolid* solid )
{
  G4String error = "Trying to divide solid " + solid->GetName()
                 + " of type " + solid->GetEntityType() + " along axis ";
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
  G4Exception("G4ReplicatedSlice::ErrorInAxis()", "GeomDiv0002",
              FatalException, error);
}

// The next methods are for specialised repeated volumes 
//     (replicas, parameterised vol.) which are completely regular.
// Currently this is not applicable to divisions  ( J.A. Nov 2005 )
// ----------------------------------------------------------------------
// IsRegularRepeatedStructure()
//
G4bool G4ReplicatedSlice::IsRegularStructure() const
{
  return false;
}           

// ----------------------------------------------------------------------
// IsRegularRepeatedStructure()
//
G4int G4ReplicatedSlice::GetRegularStructureId() const
{
  return 0;  
}           
