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
// G4Line.cc
//
// ----------------------------------------------------------------------

#include "G4Line.hh"

G4Line::G4Line (){}
G4Line::~G4Line (){}

G4Line::G4Line(const G4Line& right)
  : G4Curve(), pnt(right.pnt), dir(right.dir),
    invDir(right.invDir), v(right.v)
{
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
}

G4Line& G4Line::operator=(const G4Line& right)
{
  if (&right == this) return *this;
  
  pnt       = right.pnt;
  dir       = right.dir;
  invDir    = right.invDir;
  v         = right.v;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
  
  return *this;
}

G4Curve* G4Line::Project(const G4Transform3D& tr)
{
  G4Vector3D newDir= tr*dir;
  
  if (std::abs(newDir.x())+std::abs(newDir.y()) < kCarTolerance){
  
     newDir.setX(kCarTolerance);
     newDir.setY(kCarTolerance);
  };
  
  G4Point3D newPnt= tr*pnt;
  newDir.setZ(0);
  newPnt.setZ(0);
  
  G4Line* r= new G4Line();

  // L. Broglia : terrible mistake !!!!
  //r->Init(newDir, newPnt);
  r->Init(newPnt, newDir);

  r->SetBounds(GetPStart(), GetPEnd());
  
  return r;
}

////////////////////////////////////////////////////////////////////////////

G4bool G4Line::Tangent(G4CurvePoint&, G4Vector3D& vec)
{
  if(GetSameSense())
    vec = -dir;
  else
    vec = dir;

  return true;
}
