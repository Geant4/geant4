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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VFacet.cc 92024 2015-08-13 14:16:00Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012, M Gayer, CERN, - Reviewed optimized implementation.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4VFacet.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "G4GeometryTolerance.hh"

using namespace std;

const G4double G4VFacet::dirTolerance = 1.0E-14;
const G4double G4VFacet::kCarTolerance =
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

///////////////////////////////////////////////////////////////////////////////
//
G4VFacet::~G4VFacet()
{
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4VFacet::operator== (const G4VFacet &right) const
{
  G4double tolerance = kCarTolerance*kCarTolerance/4.0;

  if (GetNumberOfVertices() != right.GetNumberOfVertices())
    return false;
  else if ((GetCircumcentre()-right.GetCircumcentre()).mag2() > tolerance)
    return false;
  else if (std::fabs((right.GetSurfaceNormal()).dot(GetSurfaceNormal())) < 0.9999999999)
    return false;

  G4bool coincident  = true;
  G4int i = 0;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    coincident = false;
    G4int j   = 0; 
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      coincident = (GetVertex(i)-right.GetVertex(j)).mag2() < tolerance;
    } while (!coincident && ++j < GetNumberOfVertices());
  } while (coincident && ++i < GetNumberOfVertices());

  return coincident;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4VFacet::ApplyTranslation(const G4ThreeVector v)
{
  G4int n = GetNumberOfVertices();
  for (G4int i = 0; i < n; ++i)
  {
    SetVertex(i, GetVertex(i) + v);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream &G4VFacet::StreamInfo(std::ostream &os) const
{
  os << G4endl;
  os << "*********************************************************************"
     << G4endl;
  os << "FACET TYPE       = " << GetEntityType() << G4endl;
  os << "ABSOLUTE VECTORS = " << G4endl;
  G4int n = GetNumberOfVertices();
  for (G4int i = 0; i < n; ++i)
    os << "P[" << i << "]      = " << GetVertex(i) << G4endl;
  os << "*********************************************************************"
     << G4endl;

  return os;
}

G4bool G4VFacet::IsInside (const G4ThreeVector &p) const
{
  G4ThreeVector d =  p - GetVertex(0);
  G4double displacement = d.dot(GetSurfaceNormal());
  return displacement <= 0.0;
}
