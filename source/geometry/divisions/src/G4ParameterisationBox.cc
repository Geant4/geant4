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
// $Id: G4ParameterisationBox.cc 68040 2013-03-13 14:19:04Z gcosmo $
//
// class G4ParameterisationBox Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
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
    std::ostringstream message;
    message << "Only axes along X are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationBoxX::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
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

  G4double pDx = fwidth/2. - fhgap;
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
    std::ostringstream message;
    message << "Only axes along Y are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationBoxY::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
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
  G4double pDy = fwidth/2. - fhgap;
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
    std::ostringstream message;
    message << "Only axes along Z are allowed !  Axis: " << faxis;
    G4Exception("G4ParameterisationBoxZ::ComputeTransformation()",
                "GeomDiv0002", FatalException, message);
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
  G4double pDz = fwidth/2. - fhgap;

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

