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
// $Id: G4Curve.cc,v 1.6 2001-07-11 09:59:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Curve.cc
//
// ----------------------------------------------------------------------

#include "G4Curve.hh"

G4Curve::G4Curve()
 : bBox(G4BoundingBox3D::space), bounded(false), sameSense(true)
{
}

G4Curve::~G4Curve()
{
}

G4Curve::G4Curve(const G4Curve& c)
 : start(c.start), end(c.end), pStart(c.pStart), pEnd(c.pEnd),
   pRange(c.pRange), bounded(c.bounded), sameSense(c.sameSense)
{
}

G4Curve& G4Curve::operator=(const G4Curve& c)
{
  if (&c == this) return *this;
  start     = c.start;
  end       = c.end;
  pStart    = c.pStart;
  pEnd      = c.pEnd;
  pRange    = c.pRange;
  bounded   = c.bounded;
  sameSense = c.sameSense;

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
