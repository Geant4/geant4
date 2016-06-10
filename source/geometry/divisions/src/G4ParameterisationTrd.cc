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
// $Id: G4ParameterisationTrd.cc 68040 2013-03-13 14:19:04Z gcosmo $
//
// class G4ParameterisationTrd Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
// --------------------------------------------------------------------

#include "G4ParameterisationTrd.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ReflectedSolid.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"

//--------------------------------------------------------------------------
G4VParameterisationTrd::
G4VParameterisationTrd( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid ),
     bDivInTrap(false)
{
  G4Trd* msol = (G4Trd*)(msolid);
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Trd*)(mConstituentSolid);
  
    // Create a new solid with inversed parameters
    G4Trd* newSolid
      = new G4Trd(msol->GetName(),
                  msol->GetXHalfLength2(), msol->GetXHalfLength1(),
                  msol->GetYHalfLength2(), msol->GetYHalfLength1(),
                  msol->GetZHalfLength());
    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//------------------------------------------------------------------------
G4VParameterisationTrd::~G4VParameterisationTrd()
{
}

//------------------------------------------------------------------------
G4ParameterisationTrdX::
G4ParameterisationTrdX( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationTrd( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionTrdX" );

  G4Trd* msol = (G4Trd*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( msol->GetXHalfLength1()+msol->GetXHalfLength2(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( msol->GetXHalfLength1()+msol->GetXHalfLength2(),
                             nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTrdX - ## divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif

  G4double mpDx1 = msol->GetXHalfLength1();
  G4double mpDx2 = msol->GetXHalfLength2();
  if( std::fabs(mpDx1 - mpDx2) > kCarTolerance )
  {
    bDivInTrap = true;
  }
}

//------------------------------------------------------------------------
G4ParameterisationTrdX::~G4ParameterisationTrdX()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationTrdX::GetMaxParameter() const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  return (msol->GetXHalfLength1()+msol->GetXHalfLength2());
}

//------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeTransformation( const G4int copyNo,
                       G4VPhysicalVolume *physVol ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdx = ( msol->GetXHalfLength1() + msol->GetXHalfLength2() ) / 2.;

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi;
  if( !bDivInTrap )
  {
    posi = -mdx + foffset + (copyNo+0.5)*fwidth;
  }
  else
  {
    G4double aveHL = (msol->GetXHalfLength1()+msol->GetXHalfLength2())/2.;
    posi = - aveHL + foffset + (copyNo+0.5)*aveHL/fnDiv*2;
  }
  if( faxis == kXAxis )
  {
    origin.setX( posi ); 
  }
  else
  { 
    std::ostringstream message;
    message << "Only axes along X are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationTrdX::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
  }

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8)
           << " G4ParameterisationTrdX::ComputeTransformation() "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }
#endif

   //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeDimensions( G4Trd& trd, const G4int, const G4VPhysicalVolume* ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double pDy1 = msol->GetYHalfLength1();
  G4double pDy2 = msol->GetYHalfLength2();
  G4double pDz = msol->GetZHalfLength();
  G4double pDx = fwidth/2. - fhgap;
  
  trd.SetAllParameters ( pDx, pDx, pDy1, pDy2, pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTrdX::ComputeDimensions():"
           << copyNo << G4endl;
    trd.DumpInfo();
  }
#endif
}

G4VSolid* G4ParameterisationTrdX::
ComputeSolid(const G4int i, G4VPhysicalVolume * pv)
{
  if( bDivInTrap ) 
  {
    return G4VDivisionParameterisation::ComputeSolid(i, pv);
  } 
  else 
  {
    return fmotherSolid;
  }
}


//--------------------------------------------------------------------------
void
G4ParameterisationTrdX::ComputeDimensions( G4Trap& trap, const G4int copyNo,
                                           const G4VPhysicalVolume* ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  G4double pDy1 = msol->GetYHalfLength1();
  G4double pDy2 = msol->GetYHalfLength2();
  G4double pDz = msol->GetZHalfLength();
  G4double pDx1 = msol->GetXHalfLength1()/fnDiv; 
  G4double pDx2 = msol->GetXHalfLength2()/fnDiv; 

  G4double cxy1 = -msol->GetXHalfLength1() + foffset
                + (copyNo+0.5)*pDx1*2;// centre of the side at y=-pDy1
  G4double cxy2 = -msol->GetXHalfLength2() + foffset
                + (copyNo+0.5)*pDx2*2;// centre of the side at y=+pDy1
  G4double alp = std::atan( (cxy2-cxy1)/pDz );
  
  trap.SetAllParameters ( pDz,
			  0.,
			  0.,
			  pDy1,
			  pDx1,
			  pDx2,
			  alp,
			  pDy2,
			  pDx1 - fhgap,
			  pDx2 - fhgap * pDx2/pDx1,
			  alp);

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTrdX::ComputeDimensions():"
           << copyNo << G4endl;
    trap.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
void G4ParameterisationTrdX::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();
/*
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double mpDx1 = msol->GetXHalfLength1();
  G4double mpDx2 = msol->GetXHalfLength2();
  bDivInTrap = false;

  if( std::fabs(mpDx1 - mpDx2) > kCarTolerance )
  {
    std::ostringstream message;
    message << "Invalid solid specification. NOT supported." << G4endl
            << "Making a division of a TRD along axis X," << G4endl
            << "while the X half lengths are not equal," << G4endl
            << "is not (yet) supported. It will result" << G4endl
            << "in non-equal division solids.";
    G4Exception("G4ParameterisationTrdX::CheckParametersValidity()",
                "GeomDiv0001", FatalException, message);
  }
*/
}

//--------------------------------------------------------------------------
void G4ParameterisationTrdX::
ComputeTrapParams()
{
}

//--------------------------------------------------------------------------
G4ParameterisationTrdY::
G4ParameterisationTrdY( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType)
  : G4VParameterisationTrd( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionTrdY" );

  G4Trd* msol = (G4Trd*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*msol->GetYHalfLength1(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*msol->GetYHalfLength1(),
                             nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
     G4cout << " G4ParameterisationTrdY no divisions " << fnDiv
            << " = " << nDiv << G4endl
            << " Offset " << foffset << " = " << offset << G4endl
            << " width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationTrdY::~G4ParameterisationTrdY()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationTrdY::GetMaxParameter() const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  return 2*msol->GetYHalfLength1();
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdY::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdy = msol->GetYHalfLength1();

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdy + foffset + (copyNo+0.5)*fwidth;

  if( faxis == kYAxis )
  {
    origin.setY( posi ); 
  }
  else
  { 
    std::ostringstream message;
    message << "Only axes along Y are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationTrdY::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
  }

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8)
           << " G4ParameterisationTrdY::ComputeTransformation " << copyNo
           << " pos " << origin << " rot mat " << " axis " << faxis << G4endl;
  }
#endif

   //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdY::
ComputeDimensions(G4Trd& trd, const G4int, const G4VPhysicalVolume*) const
{
  //---- The division along Y of a Trd will result a Trd, only 
  //--- if Y at -Z and +Z are equal, else use the G4Trap version
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  
  G4double pDx1 = msol->GetXHalfLength1();
  G4double pDx2 = msol->GetXHalfLength2();
  G4double pDz = msol->GetZHalfLength();
  G4double pDy = fwidth/2. - fhgap;
 
  trd.SetAllParameters ( pDx1, pDx2, pDy, pDy, pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTrdY::ComputeDimensions():" << G4endl;
    trd.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
void G4ParameterisationTrdY::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double mpDy1 = msol->GetYHalfLength1();
  G4double mpDy2 = msol->GetYHalfLength2();

  if( std::fabs(mpDy1 - mpDy2) > kCarTolerance )
  {
    std::ostringstream message;
    message << "Invalid solid specification. NOT supported." << G4endl
            << "Making a division of a TRD along axis Y while" << G4endl
            << "the Y half lengths are not equal is not (yet)" << G4endl
            << "supported. It will result in non-equal" << G4endl
            << "division solids.";
    G4Exception("G4ParameterisationTrdY::CheckParametersValidity()",
                "GeomDiv0001", FatalException, message);
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTrdZ::
G4ParameterisationTrdZ( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType )
  : G4VParameterisationTrd( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivTrdZ" );

  G4Trd* msol = (G4Trd*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*msol->GetZHalfLength(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*msol->GetZHalfLength(),
                             nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTrdZ no divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationTrdZ::~G4ParameterisationTrdZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationTrdZ::GetMaxParameter() const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  return 2*msol->GetZHalfLength();
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdZ::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdz = msol->GetZHalfLength();

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdz + OffsetZ() + (copyNo+0.5)*fwidth;
  if( faxis == kZAxis )
  {
    origin.setZ( posi ); 
  }
  else
  { 
    std::ostringstream message;
    message << "Only axes along Z are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationTrdZ::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdZ: "
           << copyNo << G4endl
           << " Position: " << origin << " - Offset: " << foffset 
           << " - Width: " << fwidth << " Axis " << faxis << G4endl;
  }
#endif

   //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdZ::
ComputeDimensions(G4Trd& trd, const G4int copyNo,
                  const G4VPhysicalVolume*) const
{
  //---- The division along Z of a Trd will result a Trd
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double pDx1 = msol->GetXHalfLength1();
  G4double DDx = (msol->GetXHalfLength2() - msol->GetXHalfLength1() );
  G4double pDy1 = msol->GetYHalfLength1();
  G4double DDy = (msol->GetYHalfLength2() - msol->GetYHalfLength1() );
  G4double pDz = fwidth/2. - fhgap;
  G4double zLength = 2*msol->GetZHalfLength();
 
  trd.SetAllParameters( pDx1+DDx*(OffsetZ()+copyNo*fwidth+fhgap)/zLength,
                        pDx1+DDx*(OffsetZ()+(copyNo+1)*fwidth-fhgap)/zLength, 
                        pDy1+DDy*(OffsetZ()+copyNo*fwidth+fhgap)/zLength,
                        pDy1+DDy*(OffsetZ()+(copyNo+1)*fwidth-fhgap)/zLength,
                        pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTrdZ::ComputeDimensions()"
           << " - Mother TRD " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised TRD: "
           << copyNo << G4endl;
    trd.DumpInfo();
  }
#endif
}
