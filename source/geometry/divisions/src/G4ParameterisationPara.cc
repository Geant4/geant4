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
// $Id: G4ParameterisationPara.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// class G4ParameterisationPara Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
// --------------------------------------------------------------------

#include "G4ParameterisationPara.hh"

#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ReflectedSolid.hh"
#include "G4Para.hh"

//--------------------------------------------------------------------------
G4VParameterisationPara::
G4VParameterisationPara( EAxis axis, G4int nDiv, G4double width,
                         G4double offset, G4VSolid* msolid,
                         DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  G4Para* msol = (G4Para*)(msolid);
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Para*)(mConstituentSolid);
    fmotherSolid = msol;

    // Create a new solid with inversed parameters
    G4Para* newSolid
      = new G4Para(msol->GetName(),
                   msol->GetXHalfLength(), 
                   msol->GetYHalfLength(),
                   msol->GetZHalfLength(),
                   std::atan(msol->GetTanAlpha()),
                   pi - msol->GetSymAxis().theta(),
                   msol->GetSymAxis().phi());
   
    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//------------------------------------------------------------------------
G4VParameterisationPara::~G4VParameterisationPara()
{
}

//------------------------------------------------------------------------
G4ParameterisationParaX::
G4ParameterisationParaX( EAxis axis, G4int nDiv,
                         G4double width, G4double offset,
                         G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPara( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionParaX" );

  G4Para* mpara = (G4Para*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mpara->GetXHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
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
  G4double pDx = fwidth/2. - fhgap;
  G4double pDy = msol->GetYHalfLength();
  G4double pDz = msol->GetZHalfLength();
  G4double pAlpha = std::atan(msol->GetTanAlpha());
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
  :  G4VParameterisationPara( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionParaY" );

  G4Para* mpara = (G4Para*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mpara->GetYHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
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
  G4double pDy = fwidth/2. - fhgap;
  G4double pDz = msol->GetZHalfLength();
  G4double pAlpha = std::atan(msol->GetTanAlpha());
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
  :  G4VParameterisationPara( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionParaZ" );

  G4Para* mpara = (G4Para*)(fmotherSolid);
  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*mpara->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
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
  G4double posi = -mdz + OffsetZ() + (copyNo+0.5)*fwidth;
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
  G4double pDz = fwidth/2. - fhgap;
  G4double pAlpha = std::atan(msol->GetTanAlpha());
  G4double pTheta = msol->GetSymAxis().theta();
  G4double pPhi = msol->GetSymAxis().phi();
 
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
