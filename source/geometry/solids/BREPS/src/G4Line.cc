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
// $Id: G4Line.cc,v 1.9 2004/12/02 09:31:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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
