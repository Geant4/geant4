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
// $Id: G4ParameterisationTrd.cc,v 1.8 2003/11/19 11:51:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// class G4ParameterisationTrd Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "G4ParameterisationTrd.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"

//------------------------------------------------------------------------
G4ParameterisationTrdX::
G4ParameterisationTrdX( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionTrdX" );

  if( divType == DivWIDTH )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fnDiv = CalculateNDiv( 2*mtrd->GetXHalfLength1(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fwidth = CalculateWidth( 2*mtrd->GetXHalfLength1(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationTrdX - ## divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationTrdX::~G4ParameterisationTrdX()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationTrdX::GetMaxParameter() const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  return 2*msol->GetXHalfLength1();
}

//------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeTransformation( const G4int copyNo,
                       G4VPhysicalVolume *physVol ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdx = msol->GetXHalfLength1();

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdx + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kXAxis )
  {
    origin.setX( posi ); 
  }
  else
  { 
    G4cerr << "ERROR - G4ParameterisationTrdX::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationTrdX::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along X are allowed !");
  }

#ifdef G4DIVDEBUG
  if( verbose >= -2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdX: "
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
  G4double pDx = fwidth/2.;
 
  trd.SetAllParameters ( pDx, pDx, pDy1, pDy2, pDz );

#ifdef G4DIVDEBUG
  if( verbose >= -2 )
  {
     G4cout << " G4ParameterisationTrdX::ComputeDimensions():" << G4endl;
     trd.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
void G4ParameterisationTrdX::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double mpDx1 = msol->GetXHalfLength1();
  G4double mpDx2 = msol->GetXHalfLength2();

  if( fabs(mpDx1 - mpDx2) > kCarTolerance )
  {
    G4cerr << "ERROR - G4ParameterisationTrdX::CheckParametersValidity()"
           << G4endl
           << "        Making a division of a TRD along axis X," << G4endl
           << "        while the X half lengths are not equal," << G4endl
           << "        is not (yet) supported. It will result" << G4endl
           << "        in non-equal division solids." << G4endl;
    G4Exception("G4ParameterisationTrdX::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                "Invalid solid specification. NOT supported.");
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTrdY::
G4ParameterisationTrdY( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType)
  : G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionTrdY" );

  if( divType == DivWIDTH )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fnDiv = CalculateNDiv( 2*mtrd->GetYHalfLength1(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fwidth = CalculateWidth( 2*mtrd->GetYHalfLength1(),
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
    G4Exception("G4ParameterisationTrdY::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along Y are allowed !");
  }

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdY " << copyNo
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
  G4double pDy = fwidth/2.;
 
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

  if( fabs(mpDy1 - mpDy2) > kCarTolerance )
  {
    G4cerr << "ERROR - G4ParameterisationTrdY::CheckParametersValidity()"
           << G4endl
           << "        Making a division of a TRD along axis Y while" << G4endl
           << "        the Y half lengths are not equal is not (yet)" << G4endl
           << "        supported. It will result in non-equal" << G4endl
           << "        division solids." << G4endl;
    G4Exception("G4ParameterisationTrdY::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                "Invalid solid specification. NOT supported.");
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTrdZ::
G4ParameterisationTrdZ( EAxis axis, G4int nDiv,
                        G4double width, G4double offset,
                        G4VSolid* msolid, DivisionType divType )
  : G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{ 
  CheckParametersValidity();
  SetType( "DivTrdZ" );

  if( divType == DivWIDTH )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fnDiv = CalculateNDiv( 2*mtrd->GetZHalfLength(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fwidth = CalculateWidth( 2*mtrd->GetZHalfLength(),
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
  G4double posi = -mdz + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kZAxis )
  {
    origin.setZ( posi ); 
  }
  else
  { 
    G4Exception("G4ParameterisationTrdZ::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along Z are allowed !");
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
  G4double pDz = fwidth/2.;
  G4double zLength = 2*msol->GetZHalfLength();
 
  trd.SetAllParameters ( pDx1+DDx*(foffset+copyNo*fwidth)/zLength,
                         pDx1+DDx*(foffset+(copyNo+1)*fwidth)/zLength, 
                         pDy1+DDy*(foffset+copyNo*fwidth)/zLength,
                         pDy1+DDy*(foffset+(copyNo+1)*fwidth)/zLength, pDz );

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
