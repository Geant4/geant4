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
// $Id: G4VDivisionParameterisation.cc,v 1.3 2003-10-21 09:04:28 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
                             G4VSolid* motherSolid )
  : faxis(axis), fnDiv( nDiv), fwidth(step),
    foffset(offset), fmotherSolid( motherSolid ) 
{
  G4int verbose = 0;
  if (verbose >= 1)
  {
    G4cout << " G4VDivisionParameterisation  no divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " offset " << foffset << " = " << offset << G4endl
           << " step " << fwidth << " = " << step << G4endl;
  }
  theVoluFirstCopyNo = 1;

  CheckAxisIsValid();

}

//--------------------------------------------------------------------------
G4VDivisionParameterisation::~G4VDivisionParameterisation()
{
}

//--------------------------------------------------------------------------
void
G4VDivisionParameterisation::
CheckMotherSolid( const G4VSolid* motherSolid,
                  const G4String& solidType ) const
{
  if( motherSolid->GetEntityType() != solidType )
  {
    G4cerr << "ERROR - G4VDivisionParameterisation::CheckMotherSolid()"
           << G4endl
           << "        Incorrect solid type in call to G4PVDivision !"
           << G4endl
           << "        It is: " << motherSolid->GetEntityType()
           << ", while it should be: " << solidType << "." << G4endl;
    motherSolid->DumpInfo();
    G4Exception("Incorrect solid type for division!");
  }
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
  G4cout << " Motherdim: " <<  motherDim << ", Offset: " << offset
         << ", Width: " << width << ", Number of divisions: "
         << ( motherDim - offset ) / width << G4endl;

  return G4int( ( motherDim - offset ) / width );
}

//--------------------------------------------------------------------------
G4double
G4VDivisionParameterisation::
CalculateWidth( G4double motherDim, G4int nDiv, G4double offset ) const
{ 
  G4cout << " CalculateWidth: " << ( motherDim - offset ) / nDiv
	 << ", Motherdim: " << motherDim << ", Offset: " << offset
	 << ", Number of divisions: " << nDiv << G4endl;
  return ( motherDim - offset ) / nDiv;
}

