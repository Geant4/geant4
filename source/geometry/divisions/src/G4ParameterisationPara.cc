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
// $Id: G4ParameterisationPara.cc,v 1.6 2003/11/19 11:51:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// class G4ParameterisationPara Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "G4ParameterisationPara.hh"

#include <iomanip>

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Para.hh"

//------------------------------------------------------------------------
G4ParameterisationParaX::
G4ParameterisationParaX( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionParaX" );

  if( divType == DivWIDTH )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mpara->GetXHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mpara->GetXHalfLength(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaX - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4double G4ParameterisationParaX::GetMaxParameter() const
{
  G4Para* msol = (G4Para*)(fmotherSolid);
  return 2*msol->GetXHalfLength();
}

//------------------------------------------------------------------------
G4ParameterisationParaX::~G4ParameterisationParaX()
{
}

//------------------------------------------------------------------------
void
G4ParameterisationParaX::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
  G4double mdx = msol->GetXHalfLength( );

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdx + foffset+(copyNo+0.5)*fwidth;
  origin.setX( posi ); 
  
#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaX "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }
#endif

  //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationParaX::
ComputeDimensions(G4Para& para, const G4int,
                  const G4VPhysicalVolume*) const
{
  //---- The division along X of a Para will result a Para
  G4Para* msol = (G4Para*)(fmotherSolid);

  //---- Get
  G4double pDx = fwidth/2.;
  G4double pDy = msol->GetYHalfLength();
  G4double pDz = msol->GetZHalfLength();
  G4double pAlpha = atan(msol->GetTanAlpha());
  G4double pTheta = msol->GetSymAxis().theta();
  G4double pPhi = msol->GetSymAxis().phi();
 
  para.SetAllParameters ( pDx, pDy, pDz, pAlpha, pTheta, pPhi );

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaX::ComputeDimensions()"
           << " - Mother PARA " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised PARA: " << G4endl;
    para.DumpInfo();
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationParaY::
G4ParameterisationParaY( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionParaY" );

  if( divType == DivWIDTH )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mpara->GetYHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mpara->GetYHalfLength(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaY - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationParaY::~G4ParameterisationParaY()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationParaY::GetMaxParameter() const
{
  G4Para* msol = (G4Para*)(fmotherSolid);
  return 2*msol->GetYHalfLength();
}

//------------------------------------------------------------------------
void
G4ParameterisationParaY::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
  G4double mdy = msol->GetYHalfLength( );

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posiY = -mdy + foffset+(copyNo+0.5)*fwidth;
  origin.setY( posiY );
  G4double posiX = posiY * msol->GetTanAlpha();
  origin.setX( posiX );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaY "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }
#endif

  //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationParaY::
ComputeDimensions(G4Para& para, const G4int,
                  const G4VPhysicalVolume*) const
{
  //---- The division along Y of a Para will result a Para
  G4Para* msol = (G4Para*)(fmotherSolid);

  //---- Get
  G4double pDx = msol->GetXHalfLength();
  G4double pDy = fwidth/2.;
  G4double pDz = msol->GetZHalfLength();
  G4double pAlpha = atan(msol->GetTanAlpha());
  G4double pTheta = msol->GetSymAxis().theta();
  G4double pPhi = msol->GetSymAxis().phi();
 
  para.SetAllParameters ( pDx, pDy, pDz, pAlpha, pTheta, pPhi );

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationParaY::ComputeDimensions()"
           << " - Mother PARA " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised PARA: " << G4endl;
    para.DumpInfo();
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationParaZ::
G4ParameterisationParaZ( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionParaZ" );

  if( divType == DivWIDTH )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mpara->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mpara = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mpara->GetZHalfLength(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationParaZ - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationParaZ::~G4ParameterisationParaZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationParaZ::GetMaxParameter() const
{
  G4Para* msol = (G4Para*)(fmotherSolid);
  return 2*msol->GetZHalfLength();
}

//------------------------------------------------------------------------
void
G4ParameterisationParaZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
  G4double mdz = msol->GetZHalfLength( );

  //----- translation 
  G4double posi = -mdz + foffset + (copyNo+0.5)*fwidth;
  G4ThreeVector symAxis = msol->GetSymAxis();
  G4ThreeVector origin( symAxis * posi / symAxis.z() ); 
  
#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }
#endif

  //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationParaZ::
ComputeDimensions(G4Para& para, const G4int,
                  const G4VPhysicalVolume*) const
{
  //---- The division along Z of a Para will result a Para
  G4Para* msol = (G4Para*)(fmotherSolid);

  //---- Get
  G4double pDx = msol->GetXHalfLength();
  G4double pDy = msol->GetYHalfLength();
  G4double pAlpha = atan(msol->GetTanAlpha());
  G4double pTheta = msol->GetSymAxis().theta();
  G4double pPhi = msol->GetSymAxis().phi();
  G4double pDz = fwidth/2.;
 
  para.SetAllParameters ( pDx, pDy, pDz, pAlpha, pTheta, pPhi );

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationParaZ::ComputeDimensions()"
           << " - Mother PARA " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised PARA: " << G4endl;
    para.DumpInfo();
  }
#endif
}
