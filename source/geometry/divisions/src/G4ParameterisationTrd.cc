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
// $Id: G4ParameterisationTrd.cc,v 1.1 2003-06-16 15:11:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionTrdX" );

  if( divType == DivWIDTH )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fnDiv = CalculateNDiv( mtrd->GetXHalfLength1()+mtrd->GetXHalfLength2(),
                           width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Trd* mtrd = (G4Trd*)(msolid);
    fwidth = CalculateWidth( mtrd->GetXHalfLength1()+ mtrd->GetXHalfLength2(),
                             nDiv, offset );
  }

  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationTrdX - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationTrdX::~G4ParameterisationTrdX()
{
}

//------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeTransformation( const G4int copyNo,
                       G4VPhysicalVolume *physVol ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdx = (msol->GetXHalfLength1() + msol->GetXHalfLength2())/2.;

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
    G4Exception("G4ParameterisationTrdX - Only axes along X are allowed !");
  }

  if( verbose >= -2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdX: "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }

   //----- set translation 
  physVol->SetTranslation( origin );
}

//------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeDimensions( G4Trap& trap, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  //---- The division of a Trd will result a Trap
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  //---- Get the X at -/+pDz
  //----- X step and offset are given in the middle of the Trd.
  //----- In the -Z and +Z the extension in X is different.
  //----- Calculate the proportional step and offset at -Z and +Z
  G4double mpDx1 = msol->GetXHalfLength1();
  G4double mpDx2 = msol->GetXHalfLength2();
  G4double propX1 =  mpDx1 / ((mpDx1+mpDx2)/2);
  G4double propX2 =  mpDx2 / ((mpDx1+mpDx2)/2);
/*
  G4double xZneg = msol->GetXHalfLength1() * ( -fnDiv + 2*copyNo + 1 )
                 + foffset * propX1;
  G4double xZpos = msol->GetXHalfLength2() * ( -fnDiv + 2*copyNo + 1 )
                 + foffset * propX2;
*/
  G4double pDz = msol->GetZHalfLength();
  G4double pTheta = 0.;
  G4double pPhi = 0.;

  G4double mpDy1 = msol->GetYHalfLength1();
  G4double mpDy2 = msol->GetYHalfLength2();

  // four corners of TRD are  (-X1,-Y1), (-X2,Y1), (X2,Y2), (X1,-Y2)
  // the division centres in X are px =
  //
  G4double xd = foffset + fwidth * (copyNo+0.5);
  G4double fx = xd / ( mpDx1 + mpDx2 );
  G4double pDy1 = mpDy1 + (mpDy2-mpDy1) * fx / (mpDy2+mpDy1)/2.*2.;

  // get width at lower Y
  //
  G4double pDx1 = fwidth/2. * propX1;
  G4double pDx2 = fwidth/2. * propX2;

  // point at centre lower Y
  //
  G4double pAlp1 = atan( ( mpDx2*(2*fx-1) - mpDx1*(2*fx-1) ) 
                         / 2 / ( mpDy1*(1-fx) + mpDy2*fx) );

  G4double pDy2 = pDy1;
  G4double pDx3 = pDx1;
  G4double pDx4 = pDx2;
  G4double pAlp2 = pAlp1;
 
  trap.SetAllParameters ( pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1,
                          pDy2, pDx3, pDx4, pAlp2);

  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationTrdX::ComputeDimensions(G4Trap)"
           << " - Mother TRD: " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised TRAP: "
           << copyNo << G4endl;
    trap.DumpInfo();
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdX::
ComputeDimensions( G4Trd& trd, const G4int copyNo,
                   const G4VPhysicalVolume* pv ) const
{
  //---- The division along X of a Trd will result a Trd, only 
  //--- if X at -Z and +Z are equal, else use the G4Trap version
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double mpDx1 = msol->GetXHalfLength1();
  G4double mpDx2 = msol->GetXHalfLength2();
  if( mpDx1 != mpDx2 )
  {
    G4cerr << "ERROR - G4ParameterisationTrdX::ComputeDimensions(G4Trd)"
           << G4endl
           << "        Different values of X halfwidth at -Z: " << mpDx1
           << " and +Z " << mpDx2 << G4endl
           << "        Solid should be a G4Trap instead of a G4Trd !"
           << G4endl;
    G4Exception("G4ParameterisationTrdX - Error in ComputeDimensions()");
  }

  //---- Get the X at -/+pDz
  G4double pDy1 = msol->GetYHalfLength1();
  G4double pDy2 = msol->GetYHalfLength2();
  G4double pDz = msol->GetZHalfLength();
  G4double pDx = fwidth/2.;
 
  trd.SetAllParameters ( pDx, pDx, pDy1, pDy2, pDz );

  if( verbose >= -2 )
  {
     G4cout << " G4ParameterisationTrdX::ComputeDimensions(G4Trd)" << G4endl
            << " Name: " << pv->GetName() << " - copyNo: " << copyNo << G4endl;
     trd.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTrdY::
G4ParameterisationTrdY( EAxis axis, G4int nDiv,
                        G4double step, G4double offset,
                        G4VSolid* msolid, DivisionType)
  : G4VDivisionParameterisation( axis, nDiv, step, offset, msolid )
{
  SetType( "ReplicaTrdY" );
  if( verbose >= 1 )
  {
     G4cout << " G4ParameterisationTrdY no divisions " << fnDiv
            << " = " << nDiv << G4endl
            << " Offset " << foffset << " = " << offset << G4endl
            << " Step " << fwidth << " = " << step << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationTrdY::~G4ParameterisationTrdY()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdY::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Trd* msol = (G4Trd*)(fmotherSolid );
  G4double mdy = (msol->GetYHalfLength1() + msol->GetYHalfLength2())/2.;

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdy + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kYAxis )
  {
    origin.setY( posi ); 
  }
  else
  { 
    G4Exception("G4ParameterisationTrdY - Only axes along Y are allowed !");
  }

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdY " << copyNo
           << " pos " << origin << " rot mat " << " axis " << faxis << G4endl;
  }

   //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdY::
ComputeDimensions(G4Trap& trap, const G4int copyNo,
                  const G4VPhysicalVolume*) const
{
  //---- The division of a Trd will result a Trap
  G4Trd* msol = (G4Trd*)(fmotherSolid);
  
  //---- Get the Y at -/+pDz
  //----- Y step and offset are given in the middle of the Trd. In the -Z and +Z the extension in Y is different. Calculate the proportional step and offset at -Z and +Z
  G4double mpDy1 = msol->GetYHalfLength1();
  G4double mpDy2 = msol->GetYHalfLength2();
  G4double propY1 = mpDy1 / ((mpDy1+mpDy2)*2);
  G4double propY2 = mpDy2 / ((mpDy1+mpDy2)*2);
  G4double yZneg = msol->GetYHalfLength1() * ( -fnDiv + 2*copyNo + 1 ) + foffset * propY1;
  G4double yZpos = msol->GetYHalfLength2() * ( -fnDiv + 2*copyNo + 1 ) + foffset * propY2;

  G4double pDz = msol->GetZHalfLength();
  G4double pTheta = atan((yZpos - yZneg) / pDz);
  G4double pPhi = 90*deg;
  G4double pDy1 = fwidth/2. * propY1;
  G4double pDx1 = msol->GetXHalfLength1();
  G4double pDx2 = pDx1;
  G4double pAlp1 = 0.;
  G4double pDy2 = fwidth/2. * propY2;
  G4double pDx3 = msol->GetXHalfLength2();
  G4double pDx4 = pDx3;
  G4double pAlp2 = 0.;
 
  trap.SetAllParameters ( pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1,
                          pDy2, pDx3, pDx4, pAlp2);

  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationTrdY::ComputeDimensions(G4Trap)"
           << " - Mother TRD " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised TRAP: "
           << copyNo << G4endl;
    trap.DumpInfo();
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationTrdY::
ComputeDimensions(G4Trd& trd, const G4int copyNo,
                  const G4VPhysicalVolume * pv) const
{
  //---- The division along Y of a Trd will result a Trd, only 
  //--- if Y at -Z and +Z are equal, else use the G4Trap version
  G4Trd* msol = (G4Trd*)(fmotherSolid);

  G4double mpDy1 = msol->GetYHalfLength1();
  G4double mpDy2 = msol->GetYHalfLength2();
  if( mpDy1 != mpDy2 )
  {
    G4cerr << "ERROR - G4ParameterisationTrdY::ComputeDimensions(G4Trd)"
           << G4endl
           << "        Different values of Y halfwidth at -Z: " << mpDy1
           << " and +Z " << mpDy2 << G4endl
           << "        Solid should be a G4Trap instead of a G4Trd !"
           << G4endl;
    G4Exception("G4ParameterisationTrdY - Error in ComputeDimensions()");
  }

  //---- Get the Y at -/+pDz
  G4double pDx1 = msol->GetXHalfLength1();
  G4double pDx2 = msol->GetXHalfLength2();
  G4double pDz = msol->GetZHalfLength();
  G4double pDy = fwidth/2.;
 
  trd.SetAllParameters ( pDx1, pDx2, pDy, pDy, pDz );

  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationTrdY::ComputeDimensions(G4Trd)" << G4endl
           << " Name: " << pv->GetName() << " - copyNo: " << copyNo << G4endl;
    trd.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationTrdZ::
G4ParameterisationTrdZ( EAxis axis, G4int nDiv,
                        G4double offset, G4double step,
                        G4VSolid* motherSolid, DivisionType )
  : G4VDivisionParameterisation( axis, nDiv, offset, step, motherSolid )
{ 
  SetType( "ReplicaTrdZ" );
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTrdZ no divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Step " << fwidth << " = " << step << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationTrdZ::~G4ParameterisationTrdZ()
{
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
    G4Exception("G4ParameterisationTrdZ - Only axes along Z are allowed !");
  }

  if( verbose >= 1 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationTrdZ: "
           << copyNo << G4endl
           << " Position: " << origin << " Axis " << faxis << G4endl;
  }

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

  //---- Get the X at -/+pDz
  G4double pDx1 = msol->GetXHalfLength1();
  G4double pDx2 = msol->GetXHalfLength2();
  G4double pDy1 = msol->GetYHalfLength1();
  G4double pDy2 = msol->GetYHalfLength2();
  G4double pDz = fwidth/2.;
 
  trd.SetAllParameters ( pDx1, pDx2, pDy1, pDy2, pDz );

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationTrdZ::ComputeDimensions(G4Trd)"
           << " - Mother TRD " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised TRD: "
           << copyNo << G4endl;
    trd.DumpInfo();
  }
}
