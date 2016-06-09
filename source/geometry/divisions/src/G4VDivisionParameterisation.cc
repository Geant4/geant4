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
// $Id: G4VDivisionParameterisation.cc,v 1.8 2003/11/19 11:51:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// class G4VDivisionParameterisation Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "G4VDivisionParameterisation.hh" 
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"

G4int G4VDivisionParameterisation::verbose = 5;

//--------------------------------------------------------------------------
G4VDivisionParameterisation::
G4VDivisionParameterisation( EAxis axis, G4int nDiv,
                             G4double step, G4double offset,
                             DivisionType divType, G4VSolid* motherSolid )
  : faxis(axis), fnDiv( nDiv), fwidth(step),
    foffset(offset), fDivisionType(divType), fmotherSolid( motherSolid ) 
{
#ifdef G4DIVDEBUG
  if (verbose >= 1)
  {
    G4cout << " G4VDivisionParameterisation  no divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " offset " << foffset << " = " << offset << G4endl
           << " step " << fwidth << " = " << step << G4endl;
  }
#endif

  theVoluFirstCopyNo = 1;
}

//--------------------------------------------------------------------------
G4VDivisionParameterisation::~G4VDivisionParameterisation()
{
}

//--------------------------------------------------------------------------
void
G4VDivisionParameterisation::
ChangeRotMatrix( G4VPhysicalVolume *physVol, G4double rotZ ) const
{
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateZ( rotZ );
  //----- set rotation
  //----- delete first old rotation matrix 
  G4RotationMatrix* rmold = physVol->GetRotation();
  delete rmold;
  physVol->SetRotation(rm);
}

//--------------------------------------------------------------------------
G4int
G4VDivisionParameterisation::
CalculateNDiv( G4double motherDim, G4double width, G4double offset ) const
{
#ifdef G4DIVDEBUG
  G4cout << " G4VDivisionParameterisation::CalculateNDiv: "
         << ( motherDim - offset ) / width 
         << " Motherdim: " <<  motherDim << ", Offset: " << offset
         << ", Width: " << width << G4endl;
#endif

  return G4int( ( motherDim - offset ) / width );
}

//--------------------------------------------------------------------------
G4double
G4VDivisionParameterisation::
CalculateWidth( G4double motherDim, G4int nDiv, G4double offset ) const
{ 
#ifdef G4DIVDEBUG
  G4cout << " G4VDivisionParameterisation::CalculateWidth: "
         << ( motherDim - offset ) / nDiv
         << ", Motherdim: " << motherDim << ", Offset: " << offset
         << ", Number of divisions: " << nDiv << G4endl;
#endif

  return ( motherDim - offset ) / nDiv;
}

//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckParametersValidity()
{
  G4double maxPar = GetMaxParameter();
  CheckOffset( maxPar );
  CheckNDivAndWidth( maxPar );
}

//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckOffset( G4double maxPar )
{
  if( foffset >= maxPar )
  {
    G4cerr << "ERROR - G4VDivisionParameterisation::CheckOffset()" << G4endl
           << "        Division of solid " << fmotherSolid->GetName()
           << " has too big offset = " << G4endl
           << "        " << foffset << " > " << maxPar << " !" << G4endl;
    G4Exception("G4VDivisionParameterisation::CheckOffset()",
                "IllegalConstruct", FatalException,
                "Not supported configuration.");
  }
}

//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckNDivAndWidth( G4double maxPar )
{
  if( (fDivisionType == DivNDIVandWIDTH)
      && (foffset + fwidth*fnDiv - maxPar > kCarTolerance ) )
  {
    G4cerr << "ERROR - G4VDivisionParameterisation::CheckNDivAndWidth()"
           << G4endl
           << "        Division of solid " << fmotherSolid->GetName()
           << " has too big offset + width*nDiv = " << G4endl
           << "        " << foffset + fwidth*fnDiv << " > "
           << foffset << ". Width = "
           << G4endl
           << "        " << fwidth << ". nDiv = " << fnDiv << " !"
           << G4endl;
    G4Exception("G4VDivisionParameterisation::CheckNDivAndWidth()",
                "IllegalConstruct", FatalException,
                "Not supported configuration.");
  }
}
