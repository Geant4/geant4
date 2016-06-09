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
// $Id: G4CurveRayIntersection.cc,v 1.6 2006/06/29 18:42:06 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CurveRayIntersection.cc
//
// ----------------------------------------------------------------------

#include "G4CurveRayIntersection.hh"

const G4int G4CurveRayIntersection::dFlag= 4;

G4CurveRayIntersection::G4CurveRayIntersection()
  : r(0), d(kInfinity)
{
}

G4CurveRayIntersection::G4CurveRayIntersection(G4Curve& c0, const G4Ray& r0)
{
  Init(c0, r0);
}

G4CurveRayIntersection::~G4CurveRayIntersection()
{
}

G4CurveRayIntersection::G4CurveRayIntersection(const G4CurveRayIntersection& cr)
  : G4CurvePoint(), r(cr.r), d(cr.d)
{
  c = cr.c;
  p = cr.p;
  u = cr.u;
  notComputed = cr.notComputed;
}

G4CurveRayIntersection&
G4CurveRayIntersection::operator=(const G4CurveRayIntersection& cr)
{
  if (&cr == this) return *this;
  c = cr.c;
  p = cr.p;
  u = cr.u;
  r = cr.r;
  d = cr.d;
  notComputed = cr.notComputed;

  return *this;
}
