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
// $Id: G4VDivisionParameterisation.cc,v 1.7 2003-11-18 12:15:44 arce Exp $
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
                             DivisionType divType, G4VSolid* motherSolid )
  : faxis(axis), fnDiv( nDiv), fwidth(step),
    foffset(offset), fDivisionType(divType), fmotherSolid( motherSolid ) 
{
  if (verbose >= 1)
  {
    G4cout << " G4VDivisionParameterisation  no divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " offset " << foffset << " = " << offset << G4endl
           << " step " << fwidth << " = " << step << G4endl;
  }
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
  G4cout << "CalculateNDiv " << ( motherDim - offset ) / width 
	 << " Motherdim: " <<  motherDim << ", Offset: " << offset
         << ", Width: " << width << G4endl;

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


//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckParametersValidity()
{
  G4double maxPar = GetMaxParameter();
  //  G4cout << "G4VDivisionParameterisation::CheckParametersValidity maxPar " << maxPar << G4endl;
  CheckOffset( maxPar );
  CheckNDivAndWidth( maxPar );
}


//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckOffset( G4double maxPar )
{
  if( foffset >= maxPar ) {
    char msgcha[10];
    gcvt( foffset, 10, msgcha );
    char msgcha2[10];
    gcvt( maxPar, 10, msgcha2 );
    G4String msgstr = G4String("!! Division of solid ") + G4String(fmotherSolid->GetName()) + G4String(" has too big offset = ") + G4String(msgcha) + G4String(" > ") + G4String(msgcha2);
    G4Exception("G4VDivisionParameterisation::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                msgstr);
  }
}


//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckNDivAndWidth( G4double maxPar )
{
  if( fDivisionType == DivNDIVandWIDTH && foffset + fwidth*fnDiv - maxPar > kCarTolerance ) {
    char msgcha[10];
    gcvt( foffset + fwidth*fnDiv, 10, msgcha );
    char msgcha2[10];
    gcvt( foffset, 10, msgcha2 );
    char msgcha3[10];
    gcvt( fwidth, 10, msgcha3 );
    char msgcha4[10];
    gcvt( fnDiv, 10, msgcha4 );
    char msgcha5[10];
    gcvt( maxPar, 10, msgcha5 );
    G4String msgstr = G4String("!! Division of solid ") + fmotherSolid->GetName() + G4String(" has too big offset + width*nDiv = ") + G4String(msgcha) + G4String(" > ") + G4String(msgcha5) + G4String(" offset = ") + G4String(msgcha2) + G4String(" width = ") + G4String(msgcha3) + G4String(" nDiv = ") + G4String(msgcha4);
    G4Exception("G4VDivisionParameterisation::CheckParametersValidity()",
                "IllegalConstruct", FatalException,
                msgstr);
  }
}
