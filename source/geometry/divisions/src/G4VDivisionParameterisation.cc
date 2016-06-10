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
// $Id: G4VDivisionParameterisation.cc 92625 2015-09-09 12:34:07Z gcosmo $
//
// class G4VDivisionParameterisation Implementation file
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
// --------------------------------------------------------------------

#include "G4VDivisionParameterisation.hh" 
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ReflectedSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4AutoDelete.hh"

const G4int G4VDivisionParameterisation::verbose = 5;
G4ThreadLocal G4RotationMatrix* G4VDivisionParameterisation::fRot = 0;

//--------------------------------------------------------------------------
G4VDivisionParameterisation::
G4VDivisionParameterisation( EAxis axis, G4int nDiv,
                             G4double step, G4double offset,
                             DivisionType divType, G4VSolid* motherSolid )
  : faxis(axis), fnDiv( nDiv), fwidth(step), foffset(offset),
    fDivisionType(divType), fmotherSolid( motherSolid ), fReflectedSolid(false),
    fDeleteSolid(false), theVoluFirstCopyNo(1), fhgap(0.)
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
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

//--------------------------------------------------------------------------
G4VDivisionParameterisation::~G4VDivisionParameterisation()
{
  if (fDeleteSolid) delete fmotherSolid;
}

//--------------------------------------------------------------------------
G4VSolid* 
G4VDivisionParameterisation::
ComputeSolid( const G4int i, G4VPhysicalVolume* pv )
{
  G4VSolid* solid = G4VPVParameterisation::ComputeSolid(i, pv);
  if (solid->GetEntityType() == "G4ReflectedSolid")
  {
    solid = ((G4ReflectedSolid*)solid)->GetConstituentMovedSolid();
  }
  return solid;
}      

//--------------------------------------------------------------------------
void
G4VDivisionParameterisation::
ChangeRotMatrix( G4VPhysicalVolume *physVol, G4double rotZ ) const
{
  if (!fRot)
  {
    fRot = new G4RotationMatrix();
    G4AutoDelete::Register(fRot);
  }
  fRot->rotateZ( rotZ );
  physVol->SetRotation(fRot);
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
    std::ostringstream message;
    message << "Configuration not supported." << G4endl
            << "Division of solid " << fmotherSolid->GetName()
            << " has too big offset = " << G4endl
            << "        " << foffset << " > " << maxPar << " !";
    G4Exception("G4VDivisionParameterisation::CheckOffset()",
                "GeomDiv0001", FatalException, message);
  }
}

//--------------------------------------------------------------------------
void G4VDivisionParameterisation::CheckNDivAndWidth( G4double maxPar )
{
  if( (fDivisionType == DivNDIVandWIDTH)
      && (foffset + fwidth*fnDiv - maxPar > kCarTolerance ) )
  {
    std::ostringstream message;
    message << "Configuration not supported." << G4endl
            << "Division of solid " << fmotherSolid->GetName()
           << " has too big offset + width*nDiv = " << G4endl
           << "        " << foffset + fwidth*fnDiv << " > "
           << foffset << ". Width = "
           << G4endl
           << "        " << fwidth << ". nDiv = " << fnDiv << " !";
    G4Exception("G4VDivisionParameterisation::CheckNDivAndWidth()",
                "GeomDiv0001", FatalException, message);
  }
}

//--------------------------------------------------------------------------
G4double G4VDivisionParameterisation::OffsetZ() const
{
  // take into account reflection in the offset
  G4double offset = foffset;
  if (fReflectedSolid) offset = GetMaxParameter() - fwidth*fnDiv - foffset; 

  return offset;
}  

  
