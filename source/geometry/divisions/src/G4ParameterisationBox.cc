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
// $Id: G4ParameterisationBox.cc,v 1.7 2004/05/17 07:20:40 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// class G4ParameterisationBox Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------

#include "G4ParameterisationBox.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ReflectedSolid.hh"
#include "G4Box.hh"

//--------------------------------------------------------------------------
G4VParameterisationBox::
G4VParameterisationBox( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  G4Box* msol = (G4Box*)(msolid);
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Box*)(mConstituentSolid);
    fmotherSolid = msol;
    fReflectedSolid = true;
  }    
}

//--------------------------------------------------------------------------
G4VParameterisationBox::~G4VParameterisationBox()
{
}

//--------------------------------------------------------------------------
G4ParameterisationBoxX::
G4ParameterisationBoxX( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType )
  :  G4VParameterisationBox( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionBoxX" );

  G4Box* mbox = (G4Box*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mbox->GetXHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*mbox->GetXHalfLength(), nDiv, offset );
  }
#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationBoxX - no divisions "
           << fnDiv << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationBoxX::~G4ParameterisationBoxX()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationBoxX::GetMaxParameter() const
{
  G4Box* msol = (G4Box*)(fmotherSolid);
  return 2*msol->GetXHalfLength();
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxX::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid );
  G4double mdx = msol->GetXHalfLength( );

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdx + foffset+(copyNo+0.5)*fwidth;

  if( faxis == kXAxis )
  {
    origin.setX( posi ); 
  }
  else
  { 
    G4cerr << "ERROR - G4ParameterisationBoxX::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationBoxX::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along X are allowed !");
  }
#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationBoxX: "
           << copyNo << G4endl
           << " Position " << origin << " Axis " << faxis << G4endl;
  }
#endif
  //----- set translation 
  physVol->SetTranslation( origin );
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxX::
ComputeDimensions( G4Box& box, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid);

  G4double pDx = fwidth/2.;
  G4double pDy = msol->GetYHalfLength();
  G4double pDz = msol->GetZHalfLength();

  box.SetXHalfLength( pDx );
  box.SetYHalfLength( pDy );
  box.SetZHalfLength( pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationBoxX::ComputeDimensions()" << G4endl
           << " pDx: " << pDz << G4endl;
    box.DumpInfo();
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationBoxY::
G4ParameterisationBoxY( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType)
  :  G4VParameterisationBox( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionBoxY" );

  G4Box* mbox = (G4Box*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mbox->GetYHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*mbox->GetYHalfLength(), nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationBoxY - no divisions " << fnDiv << " = "
           << nDiv << ". Offset " << foffset << " = " << offset
           << ". Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationBoxY::~G4ParameterisationBoxY()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationBoxY::GetMaxParameter() const
{
  G4Box* msol = (G4Box*)(fmotherSolid);
  return 2*msol->GetYHalfLength();
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxY::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid);
  G4double mdy = msol->GetYHalfLength();

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdy + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kYAxis )
  {
    origin.setY( posi ); 
  }
  else
  {
    G4cerr << "ERROR - G4ParameterisationBoxY::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationBoxY::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along Y are allowed !");
  }
#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationBoxY: "
           << copyNo << G4endl
           << " Position " << origin << " Axis " << faxis << G4endl;
  }
#endif
   //----- set translation 
  physVol->SetTranslation( origin );
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxY::
ComputeDimensions( G4Box& box, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid);

  G4double pDx = msol->GetXHalfLength();
  G4double pDy = fwidth/2.;
  G4double pDz = msol->GetZHalfLength();

  box.SetXHalfLength( pDx );
  box.SetYHalfLength( pDy );
  box.SetZHalfLength( pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationBoxY::ComputeDimensions()" << G4endl
           << " pDx: " << pDz << G4endl;
    box.DumpInfo();
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationBoxZ::
G4ParameterisationBoxZ( EAxis axis, G4int nDiv, G4double width,
                        G4double offset, G4VSolid* msolid,
                        DivisionType divType )
  :  G4VParameterisationBox( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionBoxZ" );

  G4Box* mbox = (G4Box*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mbox->GetZHalfLength(), width, offset );
  }
  else if ( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*mbox->GetZHalfLength(), nDiv, offset );
  }
#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationBoxZ - no divisions " << fnDiv << " = "
           << nDiv << ". Offset " << foffset << " = " << offset
           << ". Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationBoxZ::~G4ParameterisationBoxZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationBoxZ::GetMaxParameter() const
{
  G4Box* msol = (G4Box*)(fmotherSolid);
  return 2*msol->GetZHalfLength();
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid );
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
    G4cerr << "ERROR - G4ParameterisationBoxZ::ComputeTransformation()"
           << G4endl;
    G4Exception("G4ParameterisationBoxZ::ComputeTransformation()",
                "IllegalConstruct", FatalException,
                "Only axes along Z are allowed !");
  }
#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationBoxZ: "
           << copyNo << G4endl
           << " Position " << origin << " Axis " << faxis << G4endl;
  }
#endif
   //----- set translation 
  physVol->SetTranslation( origin );
}

//------------------------------------------------------------------------
void
G4ParameterisationBoxZ::
ComputeDimensions( G4Box& box, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Box* msol = (G4Box*)(fmotherSolid);

  G4double pDx = msol->GetXHalfLength();
  G4double pDy = msol->GetYHalfLength();
  G4double pDz = fwidth/2.;

  box.SetXHalfLength( pDx );
  box.SetYHalfLength( pDy );
  box.SetZHalfLength( pDz );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationBoxZ::ComputeDimensions()" << G4endl
           << " pDx: " << pDz << G4endl;
    box.DumpInfo();
  }
#endif
}

