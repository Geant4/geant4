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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Curve.cc
//
// ----------------------------------------------------------------------

#include "G4Curve.hh"
#include "G4GeometryTolerance.hh"

G4Curve::G4Curve()
 : bBox(G4BoundingBox3D::space), pStart(0.), pEnd(0.), pRange(0.),
   bounded(false), sameSense(true)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4Curve::~G4Curve()
{
}

G4Curve::G4Curve(const G4Curve& c)
 : start(c.start), end(c.end), pStart(c.pStart), pEnd(c.pEnd),
   pRange(c.pRange), bounded(c.bounded), sameSense(c.sameSense),
   kCarTolerance(c.kCarTolerance)
{
}

G4Curve& G4Curve::operator=(const G4Curve& c)
{
  if (&c == this)  { return *this; }
  start     = c.start;
  end       = c.end;
  pStart    = c.pStart;
  pEnd      = c.pEnd;
  pRange    = c.pRange;
  bounded   = c.bounded;
  sameSense = c.sameSense;
  kCarTolerance = c.kCarTolerance;

  return *this;
}

G4String G4Curve::GetEntityType() const
{
  return G4String("G4Curve");
}

const char* G4Curve::Name() const
{
  return "G4Curve";
}

void G4Curve::SetParentSrfPtr(const G4Surface*)
{
}
