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
// $Id: G4ParameterisationPolycone.cc,v 1.6 2003-11-18 12:15:43 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4ParameterisationPolycone Implementation file
//
// 26.05.03 - P.Arce Initial version
//---------------------------------------------------------------------

#include "G4ParameterisationPolycone.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//---------------------------------------------------------------------
G4ParameterisationPolyconeRho::
G4ParameterisationPolyconeRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionPolyconeRho" );

  G4Polycone* msol = (G4Polycone*)(msolid);
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

//---------------------------------------------------------------------
G4ParameterisationPolyconeRho::~G4ParameterisationPolyconeRho()
{
}

//---------------------------------------------------------------------
void G4ParameterisationPolyconeRho::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH ){
    G4cerr << "!!!WARNING:  Solid " << msol->GetName() << " G4ParameterisationPolyconeRho::CheckParametersValidity. Division along R will be done with a width different for each solid section, WIDTH will not be used " << G4endl;
  }
  if( foffset != 0. ) {
    G4cerr << "!!!WARNING:  Solid " << msol->GetName() << " G4ParameterisationPolyconeZ::CheckParametersValidity. Division along  R will be done with a width different for each solid section, OFFSET will not be used " << G4endl;
  }

}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeRho::GetMaxParameter() const{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* original_pars = msol->GetOriginalParameters();
  return original_pars->Rmax[0] - original_pars->Rmin[0];
}


//---------------------------------------------------------------------
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

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeRho "
           << copyNo << G4endl
           << " Position: " << origin/mm
           << " - Width: " << fwidth/deg
           << " - Axis: " << faxis  << G4endl;
  }
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  G4int nZplanes = origparamMother->Num_z_planes;

  G4double width = 0.;
  for( G4int ii = 0; ii < nZplanes; ii++ )
  {
    width = CalculateWidth( origparamMother->Rmax[ii]
                          - origparamMother->Rmin[ii], fnDiv, foffset );
    origparam.Rmin[ii] = origparamMother->Rmin[ii]+foffset+width*copyNo;
    origparam.Rmax[ii] = origparamMother->Rmin[ii]+foffset+width*(copyNo+1);
  }

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconeRho::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::
G4ParameterisationPolyconePhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{ 
  CheckParametersValidity();
  SetType( "DivisionPolyconePhi" );

  G4Polycone* msol = (G4Polycone*)(msolid);
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
    G4cout << " G4ParameterisationPolyconePhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset/deg << " = " << offset/deg << G4endl
           << " Width " << fwidth/deg << " = " << width/deg << G4endl;
  }
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::~G4ParameterisationPolyconePhi()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconePhi::GetMaxParameter() const{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  return msol->GetEndPhi() - msol->GetStartPhi();
}


//---------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = foffset + copyNo*fwidth;
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconePhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol, -posi );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconePhi "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  origparam.Start_angle = origparamMother->Start_angle;
  origparam.Opening_angle = fwidth;

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconePhi::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeZ::
G4ParameterisationPolyconeZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType divType)
  : G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  CheckParametersValidity();
  SetType( "DivisionPolyconeZ" );

  G4Polycone* msol = (G4Polycone*)(msolid);
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  
  if( divType == DivWIDTH )
    {
    fnDiv =
      CalculateNDiv( origparamMother->Z_values[origparamMother->Num_z_planes-1]
		     - origparamMother->Z_values[0] , width, offset );
  }
  else if( divType == DivNDIV )
    {
    fwidth =
      CalculateNDiv( origparamMother->Z_values[origparamMother->Num_z_planes-1]
		     - origparamMother->Z_values[0] , nDiv, offset );
  }
  
  if( verbose >= 1 )
    {
    G4cout << " G4ParameterisationPolyconeZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }

}

//---------------------------------------------------------------------
G4ParameterisationPolyconeZ::~G4ParameterisationPolyconeZ()
{
}


//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeZ::GetMaxParameter() const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  return origparamMother->Z_values[origparamMother->Num_z_planes-1] - origparamMother->Z_values[0];
}

//---------------------------------------------------------------------
void G4ParameterisationPolyconeZ::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH ){
    G4cerr << "!!!WARNING:  Solid " << msol->GetName() << " G4ParameterisationPolyconeZ::CheckParametersValidity. Division along Z will be done splitting in the defined z_planes, WIDTH will not be used " << G4endl;
  }

  if( foffset != 0. ) {
    G4cerr << "!!!WARNING:  Solid " << msol->GetName() << " G4ParameterisationPolyconeZ::CheckParametersValidity. Division along Z will be done splitting in the defined z_planes, OFFSET will not be used " << G4endl;
  }

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  if( origparamMother->Num_z_planes-1 != fnDiv ) { 
    char msgcha[10];
    gcvt( origparamMother->Num_z_planes-1, 10, msgcha );
    char msgcha2[10];
    gcvt( fnDiv, 10, msgcha2 );
    G4String msgstr = G4String("Division along Z will be done splitting in the defined z_planes, that is the number of division will be ") + G4String(msgcha) + G4String(" instead of ") + G4String(msgcha2); 
    G4Exception("G4ParameterisationPolyconeZ::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                msgstr);
  }

}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4Polycone* msol = (G4Polycone*)(GetMotherSolid());

  //----- set translation: along Z axis
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4double posi = (origparamMother->Z_values[copyNo]
		   + origparamMother->Z_values[copyNo+1])/2;
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconeZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
  ChangeRotMatrix( physVol );

  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  // only for mother number of planes = 2!!
  //
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );

  G4double posi = (origparamMother->Z_values[copyNo]
		   + origparamMother->Z_values[copyNo+1])/2;

  origparam.Num_z_planes = 2;
  origparam.Z_values[0] = origparamMother->Z_values[copyNo] - posi;
  origparam.Z_values[1] = origparamMother->Z_values[copyNo+1] - posi;
  origparam.Rmin[0] = origparamMother->Rmin[copyNo];
  origparam.Rmin[1] = origparamMother->Rmin[copyNo+1];
  origparam.Rmax[0] = origparamMother->Rmax[copyNo];
  origparam.Rmax[1] = origparamMother->Rmax[copyNo+1];

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconeZ::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
}
