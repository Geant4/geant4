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
// $Id: G4ParameterisationPolycone.cc,v 1.3 2003-10-30 10:19:36 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4ParameterisationPolycone Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************


#include "G4ParameterisationPolycone.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//--------------------------------------------------------------------------
G4ParameterisationPolyconeRho::
G4ParameterisationPolyconeRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{
  SetType( "DivisionPolyconeRho" );

  G4Polycone* msol = (G4Polycone*)(msolid);
  dumpPolycone( *msol );
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( origparamMother->Rmax[0]
                         - origparamMother->Rmin[0], width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( origparamMother->Rmax[0]
                           - origparamMother->Rmin[0], nDiv, offset );
  }

  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationPolyconeRho - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationPolyconeRho::~G4ParameterisationPolyconeRho()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconeRho - copyNo: " << copyNo << G4endl
           << " foffset: " << foffset
           << " - fwidth: " << fwidth << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= -2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeRho " << copyNo
           << G4endl
           << " Position: " << origin/mm
           << " - Width: " << fwidth/deg
           << " - Axis: " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4cout << "G4ParameterisationPolyconeRho::ComputeDimensions()" << G4endl
           << copyNo << G4endl;

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  G4double nZplanes = origparamMother->Num_z_planes;
  G4int ii = 0;
  G4double width = 0.;
  for( ii = 0; ii < nZplanes; ii++ )
  {
    width = CalculateWidth( origparamMother->Rmax[ii] - origparamMother->Rmin[ii],
                            fnDiv, foffset );
    origparam.Rmin[ii] = origparamMother->Rmin[ii] + foffset + width * copyNo;
    origparam.Rmax[ii] = origparamMother->Rmin[ii] + foffset + width * (copyNo+1);
  }

  *(pcone.original_parameters) = origparam;
  pcone.Reset();

  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconeRho::ComputeDimensions()" << G4endl
           << copyNo << G4endl;
    dumpPolycone( pcone );
    pcone.DumpInfo();
  }
}

//--------------------------------------------------------------------------
G4ParameterisationPolyconePhi::
G4ParameterisationPolyconePhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  SetType( "DivisionPolyconePhi" );
  G4Polycone* msol = (G4Polycone*)(msolid);
  dumpPolycone( *msol );
  G4double deltaPhi = msol->GetEndPhi() - msol->GetStartPhi();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( deltaPhi, width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( deltaPhi, nDiv, offset );
  }

  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationPolyconePhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset/deg << " = " << offset/deg << G4endl
           << " Width " << fwidth/deg << " = " << width/deg << G4endl;
  }
}

//------------------------------------------------------------------------
G4ParameterisationPolyconePhi::~G4ParameterisationPolyconePhi()
{
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
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
    G4cout << " G4ParameterisationPolyconePhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol, -posi );

  if( verbose >= -2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconePhi " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  G4double width = CalculateWidth( origparamMother->Opening_angle,
			  fnDiv, foffset );
  origparam.Start_angle = foffset + origparamMother->Start_angle + copyNo *width;
  origparam.Opening_angle = width;

  *(pcone.original_parameters) = origparam;
  pcone.Reset();

  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconePhi::ComputeDimensions()" << G4endl
           << copyNo << G4endl;
    dumpPolycone( pcone );
      pcone.DumpInfo();
  }

}

//--------------------------------------------------------------------------
G4ParameterisationPolyconeZ::
G4ParameterisationPolyconeZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType divType)
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, msolid )
{ 
  CheckAxisIsValid();
  SetType( "DivisionPolyconeZ" );
  G4Polycone* msol = (G4Polycone*)(msolid);
  dumpPolycone( *msol );
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( origparamMother->Z_values[origparamMother->Num_z_planes-1] - origparamMother->Z_values[0] , width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateNDiv( origparamMother->Z_values[origparamMother->Num_z_planes-1] - origparamMother->Z_values[0] , nDiv, offset );
  }

  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationPolyconeZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }

}

//------------------------------------------------------------------------
G4ParameterisationPolyconeZ::~G4ParameterisationPolyconeZ()
{
}

//--------------------------------------------------------------------------
void G4ParameterisationPolyconeZ::CheckAxisIsValid()
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  G4cout << "G4ParameterisationPolyconeZ::CheckAxisIsValid() " << origparamMother->Num_z_planes << G4endl;
  if( origparamMother->Num_z_planes != 1 ) { 
    G4Exception("G4ParameterisationPolyconeZ::CheckAxisIsValid()",
                "IllegalConstruct", FatalException,
		//                "Making a division of a Polycone along axis Z while there more than one Z planes is not yet supported" );
                "Making a division of a Polycone along axis Z is not (yet) supported" );
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4Polycone* msol = (G4Polycone*)(GetMotherSolid());
  dumpPolycone( *msol );
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4double ZHalfLength = (origparamMother->Z_values[origparamMother->Num_z_planes - 1] - origparamMother->Z_values[0]);
  //----- set translation: along Z axis
  G4double posi = -ZHalfLength + foffset
                  + fwidth/2 + copyNo*fwidth;
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= -2 )
  {
    G4cout << " G4ParameterisationPolyconeZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeZ " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }

}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  //only for mother number of planes = 2!!
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical* origparam( origparamMother );
  G4double zFirst = origparamMother->Z_values[0];
  G4double zDiff = origparamMother->Z_values[1] - zFirst;
  G4double aRInner = (origparamMother->Rmin[1] - origparamMother->Rmin[0] ) / zDiff;
  G4double bRInner = (origparamMother->Rmin[1] + origparamMother->Rmin[0] ) / 2;
  G4double aROuter = (origparamMother->Rmax[1] - origparamMother->Rmax[0] ) / zDiff;
  G4double bROuter = (origparamMother->Rmax[1] + origparamMother->Rmax[0] ) / 2;
  G4double xMinusZ = zFirst + foffset + fwidth*copyNo;
  G4double xPlusZ = zFirst + foffset + fwidth*(copyNo+1);

  origparam->Z_values[0] = xMinusZ;
  origparam->Z_values[1] = xPlusZ;
  origparam->Rmin[0] = aRInner * xMinusZ + bRInner;
  origparam->Rmax[1] = aROuter * xMinusZ + bROuter;
  origparam->Rmin[0] = aRInner * xPlusZ + bRInner;
  origparam->Rmax[1] = aROuter * xPlusZ + bROuter;

  *(pcone.original_parameters) = *origparam;
  pcone.Reset();

  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconeZ::ComputeDimensions()" << G4endl
           << copyNo << G4endl;
    dumpPolycone( pcone );
      pcone.DumpInfo();
  }

}


void dumpPolycone( G4Polycone& pcone )
{
  G4PolyconeHistorical* op = pcone.GetOriginalParameters();

  G4cout << " POLYCONE " << pcone.GetName() << G4endl
	 << " Start_angle " << op->Start_angle << G4endl
	 << " Opening_angle " << op->Opening_angle << G4endl
	 << " Num_z_planes " << op->Num_z_planes << G4endl;

  for( G4int ii= 0; ii < op->Num_z_planes; ii++ ){
    G4cout << ii << " Z_values " << op->Z_values[ii] 
	   << " Rmin " << op->Rmin[ii] 
	   << " Rmax " << op->Rmax[ii] << G4endl;
  }
}
