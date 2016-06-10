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
// $Id: G4ParameterisationCons.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// class G4ParameterisationCons Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
// --------------------------------------------------------------------

#include "G4ParameterisationCons.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ReflectedSolid.hh"
#include "G4Cons.hh"

//--------------------------------------------------------------------------
G4VParameterisationCons::
G4VParameterisationCons( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  G4Cons* msol = (G4Cons*)(msolid);
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Cons*)(mConstituentSolid);
  
    // Create a new solid with inversed parameters
    G4Cons* newSolid
      = new G4Cons(msol->GetName(),
                   msol->GetInnerRadiusPlusZ(), msol->GetOuterRadiusPlusZ(),
                   msol->GetInnerRadiusMinusZ(), msol->GetOuterRadiusMinusZ(),
                   msol->GetZHalfLength(),
                   msol->GetStartPhiAngle(), msol->GetDeltaPhiAngle());
    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//------------------------------------------------------------------------
G4VParameterisationCons::~G4VParameterisationCons()
{
}

//--------------------------------------------------------------------------
G4ParameterisationConsRho::
G4ParameterisationConsRho( EAxis axis, G4int nDiv,
                           G4double width, G4double offset,
                           G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationCons( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionConsRho" );

  G4Cons* msol = (G4Cons*)(fmotherSolid);
  if( msol->GetInnerRadiusPlusZ() == 0. )
  {
    std::ostringstream message;
    message << "OuterRadiusMinusZ = 0" << G4endl 
            << "Width is calculated as that of OuterRadiusMinusZ !";
    G4Exception("G4ParameterisationConsRho::G4ParameterisationConsRho()",
                "GeomDiv1001", JustWarning, message);
  } 

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( msol->GetOuterRadiusMinusZ()
                         - msol->GetInnerRadiusMinusZ(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Cons* mconsol = (G4Cons*)(msolid);
    fwidth = CalculateWidth( mconsol->GetOuterRadiusMinusZ()
                           - mconsol->GetInnerRadiusMinusZ(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationConsRho - no divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset
           << " - Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationConsRho::~G4ParameterisationConsRho()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationConsRho::GetMaxParameter() const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);
  return msol->GetOuterRadiusMinusZ() - msol->GetInnerRadiusMinusZ();
}

//--------------------------------------------------------------------------
void
G4ParameterisationConsRho::
ComputeTransformation( const G4int, G4VPhysicalVolume *physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout <<  " G4ParameterisationConsRho " << G4endl
           << " Offset: " << foffset
           << " - Width: " << fwidth << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationConsRho" << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
void
G4ParameterisationConsRho::
ComputeDimensions( G4Cons& cons, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);

  G4double pRMin1 = msol->GetInnerRadiusMinusZ() + foffset + fwidth * copyNo;
  G4double pRMax1 = msol->GetInnerRadiusMinusZ() + foffset + fwidth * (copyNo+1);
 
  //width at Z Plus
  //- G4double fwidthPlus =
  //   fwidth * ( msol->GetOuterRadiusPlusZ()/ msol->GetInnerRadiusPlusZ())
  //-         / ( msol->GetOuterRadiusMinusZ() - msol->GetInnerRadiusMinusZ());
  G4double fwidthPlus = CalculateWidth( msol->GetOuterRadiusPlusZ()
                           - msol->GetInnerRadiusPlusZ(), fnDiv, foffset );
  G4double pRMin2 = msol->GetInnerRadiusPlusZ()
                  + foffset + fwidthPlus * copyNo;
  G4double pRMax2 = msol->GetInnerRadiusPlusZ()
                  + foffset + fwidthPlus * (copyNo+1);
  G4double pDz = msol->GetZHalfLength();

  G4double d_half_gap = fhgap * pRMax2 / pRMax1;
  //- already rotated  double pSR = foffset + copyNo*fwidth;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();;

  cons.SetInnerRadiusMinusZ( pRMin1 + fhgap );
  cons.SetOuterRadiusMinusZ( pRMax1 - fhgap );
  cons.SetInnerRadiusPlusZ( pRMin2 + d_half_gap );
  cons.SetOuterRadiusPlusZ( pRMax2 - d_half_gap );
  cons.SetZHalfLength( pDz );
  cons.SetStartPhiAngle( pSPhi, false );
  cons.SetDeltaPhiAngle( pDPhi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationConsRho::ComputeDimensions()" << G4endl
           << " pRMin: " << pRMin1 << " - pRMax: " << pRMax1 << G4endl;
    if( verbose >= 4 ) cons.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationConsPhi::
G4ParameterisationConsPhi( EAxis axis, G4int nDiv,
                           G4double width, G4double offset,
                           G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationCons( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionConsPhi" );

  G4Cons* msol = (G4Cons*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( msol->GetDeltaPhiAngle(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( msol->GetDeltaPhiAngle(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationConsPhi no divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationConsPhi::~G4ParameterisationConsPhi()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationConsPhi::GetMaxParameter() const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);
  return msol->GetDeltaPhiAngle();
}

//--------------------------------------------------------------------------
void
G4ParameterisationConsPhi::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = foffset  + copyNo*fwidth;

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationConsPhi - position: " << posi/deg << G4endl
           << " Origin: " << origin << " copyNo: " << copyNo
           << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl
           << " - Axis: " << faxis << G4endl;
  }
#endif

  ChangeRotMatrix( physVol, -posi );
}

//--------------------------------------------------------------------------
void
G4ParameterisationConsPhi::
ComputeDimensions( G4Cons& cons, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);

  G4double pRMin1 = msol->GetInnerRadiusMinusZ();
  G4double pRMax1 = msol->GetOuterRadiusMinusZ();
  G4double pRMin2 = msol->GetInnerRadiusPlusZ();
  G4double pRMax2 = msol->GetOuterRadiusPlusZ();
  G4double pDz = msol->GetZHalfLength();

  //- already rotated  double pSPhi = foffset + copyNo*fwidth;
  G4double pSPhi = foffset + msol->GetStartPhiAngle() + fhgap;
  G4double pDPhi = fwidth - 2.*fhgap;

  cons.SetInnerRadiusMinusZ( pRMin1 );
  cons.SetOuterRadiusMinusZ( pRMax1 );
  cons.SetInnerRadiusPlusZ( pRMin2 );
  cons.SetOuterRadiusPlusZ( pRMax2 );
  cons.SetZHalfLength( pDz );
  cons.SetStartPhiAngle( pSPhi, false );
  cons.SetDeltaPhiAngle( pDPhi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationConsPhi::ComputeDimensions" << G4endl
           << " pSPhi: " << pSPhi << " - pDPhi: " << pDPhi << G4endl;
    if( verbose >= 4 ) cons.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationConsZ::
G4ParameterisationConsZ( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationCons( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionConsZ" );

  G4Cons* msol = (G4Cons*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*msol->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*msol->GetZHalfLength(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationConsZ: # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl
           << " - Axis: " << faxis << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationConsZ::~G4ParameterisationConsZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationConsZ::GetMaxParameter() const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);
  return 2*msol->GetZHalfLength();
}

//--------------------------------------------------------------------------
void
G4ParameterisationConsZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  //----- set translation: along Z axis
  G4Cons* motherCons = (G4Cons*)(GetMotherSolid());
  G4double posi = - motherCons->GetZHalfLength() + OffsetZ()
                  + fwidth/2 + copyNo*fwidth;
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationConsZ::ComputeTransformation()" << G4endl
           << " Origin: " << origin << " - copyNo: " << copyNo << G4endl
           << " foffset: " << foffset << " - fwidth: " << fwidth
           << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );
}


//--------------------------------------------------------------------------
void
G4ParameterisationConsZ::
ComputeDimensions( G4Cons& cons, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Cons* msol = (G4Cons*)(fmotherSolid);

  G4double mHalfLength = msol->GetZHalfLength() - fhgap;
  G4double aRInner = (msol->GetInnerRadiusPlusZ()
                   - msol->GetInnerRadiusMinusZ()) / (2*mHalfLength);
  G4double bRInner = (msol->GetInnerRadiusPlusZ()
                   + msol->GetInnerRadiusMinusZ()) / 2;
  G4double aROuter = (msol->GetOuterRadiusPlusZ()
                   - msol->GetOuterRadiusMinusZ()) / (2*mHalfLength);
  G4double bROuter = (msol->GetOuterRadiusPlusZ()
                   + msol->GetOuterRadiusMinusZ()) / 2;
  G4double xMinusZ = -mHalfLength + OffsetZ() + fwidth*copyNo + fhgap;
  G4double xPlusZ  = -mHalfLength + OffsetZ() + fwidth*(copyNo+1) - fhgap;
  cons.SetInnerRadiusMinusZ( aRInner * xMinusZ + bRInner );
  cons.SetOuterRadiusMinusZ( aROuter * xMinusZ + bROuter );
  cons.SetInnerRadiusPlusZ( aRInner * xPlusZ + bRInner );
  cons.SetOuterRadiusPlusZ( aROuter * xPlusZ + bROuter );
 
  G4double pDz = fwidth / 2. - fhgap;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();

  cons.SetZHalfLength( pDz );
  cons.SetStartPhiAngle( pSPhi, false );
  cons.SetDeltaPhiAngle( pDPhi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationConsZ::ComputeDimensions()" << G4endl
           << " pDz: " << pDz << G4endl;
    if( verbose >= 4 ) cons.DumpInfo();
  }
#endif

}
