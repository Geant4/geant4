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
// $Id: G4ParameterisationPolyhedra.cc,v 1.4 2003-11-05 17:42:28 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4ParameterisationPolyhedra Implementation file
//
// 14.10.03 - P.Arce Initial version
//---------------------------------------------------------------------

#include "G4ParameterisationPolyhedra.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraRho::
G4ParameterisationPolyhedraRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionPolyhedraRho" );

  G4Polyhedra* msol = (G4Polyhedra*)(msolid);
  G4PolyhedraHistorical* original_pars = msol->GetOriginalParameters();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( original_pars->Rmax[0]
                         - original_pars->Rmin[0], width, offset );
  }
  else if( divType == DivNDIV )
  {
    G4cout  << " original_pars " << original_pars << G4endl;
    G4cout  << " original_pars rmax size " << original_pars->Rmax << " "
            <<  original_pars->Rmin << G4endl;

    fwidth = CalculateWidth( original_pars->Rmax[0]
                           - original_pars->Rmin[0], nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyhedraRho - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationPolyhedraRho::~G4ParameterisationPolyhedraRho()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraRho::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.);

  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraRho - copyNo: " << copyNo << G4endl
           << " foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraRho "
           << copyNo << G4endl
           << " Position: " << origin
           << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraRho::
ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  G4PolyhedraHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyhedraHistorical origparam( *origparamMother );
  G4int nZplanes = origparamMother->Num_z_planes;

  G4double width = 0.;
  for( G4int ii = 0; ii < nZplanes; ii++ )
  {
    width = CalculateWidth( origparamMother->Rmax[ii]
                          - origparamMother->Rmin[ii], fnDiv, foffset );
    origparam.Rmin[ii] = origparamMother->Rmin[ii]+foffset+width*copyNo;
    origparam.Rmax[ii] = origparamMother->Rmin[ii]+foffset+width*(copyNo+1);
  }

  phedra.SetOriginalParameters(&origparam); // copy values & transfer pointers
  phedra.Reset();                           // reset to new solid parameters

  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyhedraRho::ComputeDimensions()" << G4endl
           << "-- Parametrised phedra copy-number: " << copyNo << G4endl;
    phedra.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraPhi::
G4ParameterisationPolyhedraPhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  SetType( "DivisionPolyhedraPhi" );
  G4Polyhedra* msol = (G4Polyhedra*)(msolid);
  G4double deltaPhi = msol->GetEndPhi() - msol->GetStartPhi();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( deltaPhi, width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( deltaPhi, nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyhedraPhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationPolyhedraPhi::~G4ParameterisationPolyhedraPhi()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraPhi::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = foffset  + copyNo*fwidth;
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraPhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol, -posi );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraPhi " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraPhi::
ComputeDimensions( G4Polyhedra&, const G4int,
                   const G4VPhysicalVolume* ) const
{
/*
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);
  G4PolyhedraHistorical* orig_pars = msol->GetOriginalParameters();

  G4double pRMin = msol->GetInnerRadius();
  G4double pRMax = msol->GetOuterRadius();
  G4double pDz   = msol->GetZHalfLength();

  //- already rotated  double pSPhi = foffset + copyNo*fwidth;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = fwidth;

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraPhi::ComputeDimensions()" << G4endl
           << " pSPhi: " << pSPhi << " - pDPhi: " << pDPhi << G4endl;
    tubs.DumpInfo();
  }
*/
}

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraZ::
G4ParameterisationPolyhedraZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  SetType( "DivisionPolyhedraZ" );
/*
  G4Polyhedra* msol = (G4Polyhedra*)(msolid);

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( 2*msol->GetZHalfLength(), width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( 2*msol->GetZHalfLength(), nDiv, offset );
  }

  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyhedraZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
*/
}

//------------------------------------------------------------------------
G4ParameterisationPolyhedraZ::~G4ParameterisationPolyhedraZ()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraZ::
ComputeTransformation( const G4int, G4VPhysicalVolume* ) const
{
/*
  //----- set translation: along Z axis
  G4Polyhedra* motherPolyhedra = (G4Polyhedra*)(GetMotherSolid());
  G4double posi = -motherPolyhedra->GetZHalfLength() + foffset
                  + fwidth/2 + copyNo*fwidth;
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraZ " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
*/
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraZ::
ComputeDimensions( G4Polyhedra&, const G4int,
                   const G4VPhysicalVolume* ) const
{
/*
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  G4double pRMin = msol->GetInnerRadius();
  G4double pRMax = msol->GetOuterRadius();
  //  G4double pDz = msol->GetZHalfLength() / GetNoDiv();
  G4double pDz = fwidth/2.;
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraZ::ComputeDimensions - pDz: "
           << pDz << G4endl;
    tubs.DumpInfo();
  }
*/
}
