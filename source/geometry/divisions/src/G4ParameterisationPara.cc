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
// $Id: G4ParameterisationPara.cc,v 1.2 2003-10-16 10:42:42 arce Exp $
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

  if( verbose >= -1 )
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
  G4double posi = -mdz + foffset + (copyNo+0.5)*fwidth;
  G4ThreeVector symAxis = msol->GetSymAxis();
  G4ThreeVector origin( symAxis * posi / symAxis.z() ); 
  
  if( faxis != kZAxis )
  { 
    G4cerr << "ERROR - G4ParameterisationParaZ::ComputeTransformation()"
           << G4endl
           << "        Axis is along " << faxis << " !" << G4endl;
    G4Exception("G4ParameterisationParaZ - Only axes along Z are allowed !");
  }

  if( verbose >= -2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationParaZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Axis: " << faxis << G4endl;
  }

  //----- set translation 
  physVol->SetTranslation( origin );
}

//--------------------------------------------------------------------------
void
G4ParameterisationParaZ::
ComputeDimensions(G4Para& para, const G4int copyNo,
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

  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationParaZ::ComputeDimensions(G4Para)"
           << " - Mother PARA " << G4endl;
    msol->DumpInfo();
    G4cout << " - Parameterised PARA: "
           << copyNo << G4endl;
    para.DumpInfo();
  }
}
