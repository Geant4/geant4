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
// $Id: G4ParameterisationPara.cc,v 1.1 2003-06-16 15:11:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

//--------------------------------------------------------------------------
G4ParameterisationParaX::
G4ParameterisationParaX( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionParaX" );

  if( divType == DivWIDTH )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mbox->GetXHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mbox->GetXHalfLength(), nDiv, offset );
  }
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaX - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationParaX::~G4ParameterisationParaX()
{
}

//------------------------------------------------------------------------
void
G4ParameterisationParaX::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
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
    G4cerr << "ERROR - G4ParameterisationParaX::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationParaX - Only axes along X are allowed !");
  }

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaX "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }

  //----- set translation 
  physVol->SetTranslation( origin );
}

//------------------------------------------------------------------------
G4ParameterisationParaY::
G4ParameterisationParaY( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionParaY" );

  if( divType == DivWIDTH )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mbox->GetYHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mbox->GetYHalfLength(), nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaY - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationParaY::~G4ParameterisationParaY()
{
}

//------------------------------------------------------------------------
void
G4ParameterisationParaY::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
  G4double mdy = msol->GetYHalfLength( );

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdy + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kYAxis )
  {
    origin.setY( posi ); 
  }
  else
  { 
    G4cerr << "ERROR - G4ParameterisationParaY::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationParaY - Only axes along Y are allowed !");
  }

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaY "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }

  //----- set translation 
  physVol->SetTranslation( origin );
}


//------------------------------------------------------------------------
G4ParameterisationParaZ::
G4ParameterisationParaZ( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionParaZ" );

  if( divType == DivWIDTH )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fnDiv = CalculateNDiv( 2*mbox->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4Para* mbox = (G4Para*)(msolid);
    fwidth = CalculateWidth( 2*mbox->GetZHalfLength(), nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationParaZ - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationParaZ::~G4ParameterisationParaZ()
{
}

//------------------------------------------------------------------------
void
G4ParameterisationParaZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  G4Para* msol = (G4Para*)(fmotherSolid );
  G4double mdz = msol->GetZHalfLength( );

  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  G4double posi = -mdz + foffset + (copyNo+0.5)*fwidth;
  if( faxis == kZAxis )
  {
    origin.setZ( posi ); 
  }
  else
  { 
    G4cerr << "ERROR - G4ParameterisationParaZ::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationParaZ - Only axes along Z are allowed !");
  }

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }

  //----- set translation 
  physVol->SetTranslation( origin );
}
