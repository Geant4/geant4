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
// $Id: G4ParameterisationTubs.cc,v 1.3 2003-10-30 10:19:36 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4ParameterisationTubs Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "G4ParameterisationTubs.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"

//--------------------------------------------------------------------------
G4ParameterisationTubsRho::
G4ParameterisationTubsRho( EAxis axis, G4int nDiv,
                           G4double width, G4double offset,
                           G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionTubsRho" );

  if( divType == DivWIDTH )
  {
    G4Tubs* msol = (G4Tubs*)(msolid);
    fnDiv = CalculateNDiv( msol->GetOuterRadius() - msol->GetInnerRadius(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4cout << " calculate width "<< G4endl;
    G4Tubs* msol = (G4Tubs*)(msolid);
    fwidth = CalculateWidth( msol->GetOuterRadius() - msol->GetInnerRadius(),
                             nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTubsRho - no divisions " << fnDiv << " = "
           << nDiv << G4endl
	   << " Offset " << foffset << " = " << offset << G4endl
	   << " Width " << fwidth << " = " << width << G4endl
	   << " DivType " << divType << G4endl;
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTubsRho::~G4ParameterisationTubsRho()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsRho::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsRho - copyNo: " << copyNo << G4endl
           << " Offset: " << foffset/deg
           << " - Width: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTubsRho " << copyNo
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsRho::
ComputeDimensions( G4Tubs& tubs, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Tubs* msol = (G4Tubs*)(fmotherSolid);

  G4double pRMin = msol->GetInnerRadius() + fwidth * copyNo;
  G4double pRMax = msol->GetInnerRadius() + fwidth * (copyNo+1);
  G4double pDz = msol->GetZHalfLength();
  //- already rotated  G4double pSR = foffset + copyNo*fwidth;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();;

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsRho::ComputeDimensions()" << G4endl
           << " pRMin: " << pRMin << " - pRMax: " << pRMax << G4endl;
    tubs.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTubsPhi::
G4ParameterisationTubsPhi( EAxis axis, G4int nDiv,
                           G4double width, G4double offset,
                           G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  SetType( "DivisionTubsPhi" );

  if( divType == DivWIDTH )
  {
    G4Tubs* msol = (G4Tubs*)(msolid);
    fnDiv = CalculateNDiv( msol->GetDeltaPhiAngle(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Tubs* msol = (G4Tubs*)(msolid);
    fwidth = CalculateWidth( msol->GetDeltaPhiAngle(), nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTubsPhi no divisions " << fnDiv << " = "
           << nDiv << G4endl
	   << " Offset " << foffset << " = " << offset << G4endl
	   << " Width " << fwidth << " = " << width << G4endl;
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTubsPhi::~G4ParameterisationTubsPhi()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsPhi::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = foffset  + copyNo*fwidth;
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsPhi - position: " << posi/deg << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol, -posi );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTubsPhi " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsPhi::
ComputeDimensions( G4Tubs& tubs, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Tubs* msol = (G4Tubs*)(fmotherSolid);

  G4double pRMin = msol->GetInnerRadius();
  G4double pRMax = msol->GetOuterRadius();
  G4double pDz = msol->GetZHalfLength();
  //- already rotated  G4double pSPhi = foffset + copyNo*fwidth;
  G4double pSPhi = foffset + msol->GetStartPhiAngle();
  G4double pDPhi = fwidth;

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsPhi::ComputeDimensions" << G4endl
           << " pSPhi: " << pSPhi << " - pDPhi: " << pDPhi << G4endl;
    tubs.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTubsZ::
G4ParameterisationTubsZ( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  : G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  SetType( "DivisionTubsZ" );
  if( divType == DivWIDTH )
  {
    G4Tubs* msol = (G4Tubs*)(msolid);
    fnDiv = CalculateNDiv( 2*msol->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Tubs* msol = (G4Tubs*)(msolid);
    fwidth = CalculateWidth( 2*msol->GetZHalfLength(), nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTubsZ: # divisions " << fnDiv << " = "
           << nDiv << G4endl
	   << " Offset " << foffset << " = " << offset << G4endl
	   << " Width " << fwidth << " = " << width << G4endl;
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTubsZ::~G4ParameterisationTubsZ()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsZ::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  //----- set translation: along Z axis
  G4Tubs* motherTubs = reinterpret_cast<G4Tubs*>(GetMotherSolid());
  G4double posi = - motherTubs->GetZHalfLength() + foffset
                  + fwidth/2 + copyNo*fwidth;
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsZ::ComputeTransformation()" << G4endl
           << " Position: " << posi << " - copyNo: " << copyNo << G4endl
           << " foffset " << foffset/deg << " - fwidth " << fwidth/deg
           << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTubsZ " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl; 
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationTubsZ::
ComputeDimensions( G4Tubs& tubs, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Tubs* msol = (G4Tubs*)(fmotherSolid);

  G4double pRMin = msol->GetInnerRadius();
  G4double pRMax = msol->GetOuterRadius();
  //  G4double pDz = msol->GetZHalfLength() / GetNoDiv();
  G4double pDz = fwidth/2.;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTubsZ::ComputeDimensions()" << G4endl
           << " pDz: " << pDz << G4endl;
    tubs.DumpInfo();
  }
}
