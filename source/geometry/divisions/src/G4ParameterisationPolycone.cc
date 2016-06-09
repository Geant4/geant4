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
// $Id: G4ParameterisationPolycone.cc,v 1.15 2006/06/29 18:18:44 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// class G4ParameterisationPolycone Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
//---------------------------------------------------------------------

#include "G4ParameterisationPolycone.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ReflectedSolid.hh"

//-----------------------------------------------------------------------
G4VParameterisationPolycone::
G4VParameterisationPolycone( EAxis axis, G4int nDiv, G4double width,
                             G4double offset, G4VSolid* msolid,
                             DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  G4Polycone* msol = (G4Polycone*)(msolid);
  if ((msolid->GetEntityType() != "G4ReflectedSolid") && (msol->IsGeneric()))
  {
    G4String message =
        "Sorry, generic construct for G4Polycone NOT supported.\n Solid: "
      + msol->GetName();
    G4Exception("G4VParameterisationPolycone::G4VParameterisationPolycone()",
                "NotSupported", FatalException, message);
  }
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Polycone*)(mConstituentSolid);
  
    // Get parameters
    G4int   nofZplanes = msol->GetOriginalParameters()->Num_z_planes;
    G4double* zValues  = msol->GetOriginalParameters()->Z_values;
    G4double* rminValues  = msol->GetOriginalParameters()->Rmin;
    G4double* rmaxValues  = msol->GetOriginalParameters()->Rmax;

    // Invert z values
    G4double* zValuesRefl = new double[nofZplanes];
    for (G4int i=0; i<nofZplanes; i++) zValuesRefl[i] = - zValues[i];
    
    G4Polycone* newSolid
      = new G4Polycone(msol->GetName(),
                       msol->GetStartPhi(), 
                       msol->GetEndPhi() - msol->GetStartPhi(),
                       nofZplanes, zValuesRefl, rminValues, rmaxValues);

    delete [] zValuesRefl;       

    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//---------------------------------------------------------------------
G4VParameterisationPolycone::~G4VParameterisationPolycone()
{
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeRho::
G4ParameterisationPolyconeRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionPolyconeRho" );

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
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

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationPolyconeRho - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
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

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    G4cerr << "WARNING - "
           << "G4ParameterisationPolyconeRho::CheckParametersValidity()"
           << G4endl
           << "          Solid " << msol->GetName() << G4endl
           << "          Division along R will be done with a width "
           << "different for each solid section." << G4endl
           << "          WIDTH will not be used !" << G4endl;
  }
  if( foffset != 0. )
  {
    G4cerr << "WARNING - "
           << "G4ParameterisationPolyconeRho::CheckParametersValidity()"
           << G4endl
           << "          Solid " << msol->GetName() << G4endl
           << "          Division along  R will be done with a width "
           << "different for each solid section." << G4endl
           << "          OFFSET will not be used !" << G4endl;
  }

}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeRho::GetMaxParameter() const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* original_pars = msol->GetOriginalParameters();
  return original_pars->Rmax[0] - original_pars->Rmin[0];
}


//---------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
ComputeTransformation( const G4int, G4VPhysicalVolume* physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconeRho " << G4endl
           << " foffset: " << foffset
           << " - fwidth: " << fwidth << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeRho "
           << G4endl
           << " Position: " << origin/mm
           << " - Width: " << fwidth/deg
           << " - Axis: " << faxis  << G4endl;
  }
#endif
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

#ifdef G4DIVDEBUG
  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconeRho::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::
G4ParameterisationPolyconePhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionPolyconePhi" );

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4double deltaPhi = msol->GetEndPhi() - msol->GetStartPhi();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( deltaPhi, width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( deltaPhi, nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyconePhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset/deg << " = " << offset/deg << G4endl
           << " Width " << fwidth/deg << " = " << width/deg << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::~G4ParameterisationPolyconePhi()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconePhi::GetMaxParameter() const
{
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

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconePhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol, -posi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconePhi "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeDimensions( G4Polycone& pcone, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  origparam.Start_angle = origparamMother->Start_angle;
  origparam.Opening_angle = fwidth;

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconePhi::ComputeDimensions():" << G4endl;
    pcone.DumpInfo();
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeZ::
G4ParameterisationPolyconeZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType divType)
  : G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType )
{

  CheckParametersValidity();
  SetType( "DivisionPolyconeZ" );

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
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
  
#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyconeZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
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
  return std::abs (origparamMother->Z_values[origparamMother->Num_z_planes-1]
             -origparamMother->Z_values[0]);
}

//---------------------------------------------------------------------
void G4ParameterisationPolyconeZ::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    G4cerr << "WARNING - "
           << "G4ParameterisationPolyconeZ::CheckParametersValidity()"
           << G4endl
           << "          Solid " << msol->GetName() << G4endl
           << "          Division along Z will be done splitting in the "
           << "defined z_planes." << G4endl
           << "          WIDTH will not be used !" << G4endl;
  }

  if( foffset != 0. )
  {
    G4cerr << "WARNING - "
           << "G4ParameterisationPolyconeZ::CheckParametersValidity()"
           << G4endl
           << "          Solid " << msol->GetName() << G4endl
           << "          Division along Z will be done splitting in the "
           << "defined z_planes." << G4endl
           << "          OFFSET will not be used !" << G4endl;
  }

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  if( origparamMother->Num_z_planes-1 != fnDiv )
  { 
    G4cerr << "ERROR - "
           << "G4ParameterisationPolyconeZ::CheckParametersValidity()"
           << G4endl
           << "        Division along Z will be done splitting in the defined"
           << G4endl
           << "        z_planes, i.e, the number of division would be :"
           << "        " << origparamMother->Num_z_planes-1
           << " instead of " << fnDiv << " !"
           << G4endl; 
    G4Exception("G4ParameterisationPolyconeZ::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                "Not supported configuration.");
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

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconeZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
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

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconeZ::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
#endif
}
