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
// $Id: G4CurvePoint.cc,v 1.4 2001-07-11 09:59:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
 : c(0)
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
