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
// G4CurvePoint.cc
//
// ----------------------------------------------------------------------

#include "G4CurvePoint.hh"

const G4int G4CurvePoint::pFlag= 1;
const G4int G4CurvePoint::uFlag= 2;
const G4int G4CurvePoint::allFlags= 0xFF; // lots of bits...

G4CurvePoint::G4CurvePoint()
 : c(0), u(0.), notComputed(allFlags)
{
}

G4CurvePoint::G4CurvePoint(G4Curve& c0)
{
  Init(c0);
}

G4CurvePoint::~G4CurvePoint()
{
}

G4CurvePoint::G4CurvePoint(const G4CurvePoint& cp)
  : c(cp.c), p(cp.p), u(cp.u), notComputed(cp.notComputed)
{
}

G4CurvePoint& G4CurvePoint::operator=(const G4CurvePoint& cp)
{
  if (&cp == this) return *this;
  c = cp.c;
  p = cp.p;
  u = cp.u;
  notComputed = cp.notComputed;
  
  return *this;
}
